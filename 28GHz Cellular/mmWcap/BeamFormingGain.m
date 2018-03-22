classdef BeamFormingGain < hgsetget
    properties
        nrx = 16;
        ntx = 16;
        nsample = 500;
        fc = 28;    % GHz -- used in imaging only
        
        ncluster;   % array
        Ph;         % array
        phi;        % cell
        thetA;      % cell
        thetB;      % cell
        del;        % cell
        
        
        % lamNcluster = 2.5 % 1.85; % not used
        % stdDel = 10.7*pi/180;
        meanDel = 7.8*pi/180;
        % statistics for each cluster
        stdPhi = 8.36; % dB
        %muPhi = 0.5;
        
        u;              % optimum kronecker Rx side bf vector
        %v;              % optimum kronecker Tx side bf vector
        eigProdKron;           % kronecker model gain
        gainStatDes;         % optimum kronecker bf vectors applied on channel statistics
        
        skipSecComp = true;    % make it false to compute gains wrto the second model
%         u2;             %
%         v2;             %
%         gainStatDes2;          % optimum bf vectors wrto new model applied on channel statistics
%         gainDes2;     % optimum bf vectors wrto new model applied on channel realizations
        
        gainCohDes;        % optimum bf vectors computed for every realization
        
        gainStatInt;          % optimum kronecker vectors applied on channel statistics with different random parameters
        
        minComp = false;     % make it true to compute the ones below ONLY:
        % THE ONES WE USED IN SIMULATIONS
        gainInt;         % optimum kronecker vectors applied on random channel realizations with different parameters
        gainDes;      % optimum kronecker bf vectors applied on channel realizations
        
    end
    
    methods
        function obj = BeamFormingGain(ntx,nrx)
            
            if nargin>0
                obj.ntx = ntx;
            end
            if nargin>1
                obj.nrx = nrx;
            end
            
            
            % randomly generated parameters
            
            % obj.ncluster = ceil(14*rand(obj.nsample,1));
            % obj.ncluster =ceil( -obj.lamNcluster.* log(rand(obj.nsample,1)));
            % obj.ncluster = ceil(exprnd(obj.lamNcluster,obj.nsample,1));
            obj.ncluster = 3*ones(obj.nsample,1); % constant
            obj.Ph = ones(obj.nsample,1);   % constant
            obj.del = cell(obj.nsample,1);
            obj.thetA = cell(obj.nsample,1);
            obj.thetB = cell(obj.nsample,1);
            obj.phi = cell(obj.nsample,1);
            for n = 1 : obj.nsample
                obj.del{n} = -obj.meanDel.* log(rand(obj.ncluster(n),1));
                for m = 1 :obj.ncluster(n)
                    while obj.del{n}(m)>pi
                        obj.del{n}(m) = obj.del{n}(m)-pi;
                    end
                end
                obj.thetA{n} = 2*pi*rand(obj.ncluster(n),1);
                obj.thetB{n} = 2*pi*rand(1,1);
                
                % obj.phi{n} = ones(obj.ncluster(n),1);
                %obj.phi{n} = exp(obj.stdPhi*randn(obj.ncluster(n),1));
                obj.phi{n} = 10.^(-0.1*(obj.stdPhi*randn(obj.ncluster(n),1)));
                obj.phi{n} = obj.phi{n}/sum(obj.phi{n});    % normalization
            end
            
            obj.cmpBfGain();
            %obj.plotData();
        end
        
        function cmpBfGain(obj)
            obj.gainDes = zeros(obj.nsample,1);
            obj.gainInt = zeros(obj.nsample,1);
            obj.u = cell(obj.nsample,1);
            %obj.v = cell(obj.nsample,1);
            
            if ~obj.minComp
                obj.gainStatDes = zeros(obj.nsample,1);
                obj.gainCohDes = zeros(obj.nsample,1);
                obj.gainStatInt = zeros(obj.nsample,1);
                obj.eigProdKron = zeros(obj.nsample,1);
            end
            
            if ~obj.skipSecComp
                obj.gainStatDes2 = zeros(obj.nsample,1);
                obj.gainDes2 = zeros(obj.nsample,1);
                obj.u2 = cell(obj.nsample,1);
                %obj.v2 = cell(obj.nsample,1);
            end
            
            for n = 1 : obj.nsample
                
                nsub = 20; % 20 subpath for each cluster
                m = ceil((obj.nsample-1)*rand(1));      % select a random parameter set to test for interference
                if m >= n
                    m = m+1;                            % be sure m~=n
                end
                
                
                QAk = obj.oneSideCor(obj.thetA{n}, obj.del{n}, obj.nrx);
                %QBk = obj.oneSideCor(obj.thetB{n}, obj.del{n}, obj.ntx);
                QA = zeros(obj.nrx);
                %QB = zeros(obj.ntx);
                
                for k = 1 : obj.ncluster(n)
                    QA = QA + QAk{k}*obj.phi{n}(k);
                    %QB = QB + QBk{k}*obj.phi{n}(k);
                end
                
                %RAi = (obj.Ph(n)*QA);
                %RBi = (obj.Ph(n)*QB);
                
                %RAi05 = (RAi)^0.5;
                %RBi05 = (RBi)^0.5;
                
                %[V,Diag] = eig(RAi05*(RAi05'));
                [V,Diag] = eig(QA);
                [lam1,tempind] = max(real(diag(Diag)));
                ui = V(:,tempind);
                
                %---------tx gain begin---------------
                QB = obj.oneSideCor(obj.thetB{n}, obj.del{n}(ceil(rand(1)*3)), obj.ntx);
                QB = QB{1};
                [~,DiagQB] = eig(QB);
                [txGain,~] = max(real(diag(DiagQB)));
                %---------tx gain end---------------
                
                lam2 = txGain;
                %lam2 = obj.ntx;
                %[V,Diag] = eig((RBi05')*(RBi05));
                %[lam2,tempind] = max(real(diag(Diag)));
                %vi = V(:,tempind);
                
                obj.u{n} = ui;
                %obj.v{n} = vi;
                
                
                
                if ~obj.minComp
                    gaini = lam1*lam2/obj.Ph(n);
                    obj.eigProdKron(n) = gaini;
                    
                    gainStatDesi = 0;
                    for k = 1:obj.ncluster(n)
                        gainStatDesi = gainStatDesi + real(obj.phi{n}(k)*(ui'*(QAk{k})*ui)*txGain);
                        %gainStatDesi = gainStatDesi + real(obj.phi{n}(k)*(ui'*(QAk{k})*ui)*obj.ntx);
                        %gainStatDesi = gainStatDesi + real(obj.phi{n}(k)*(ui'*(QAk{k})*ui)*(vi'*QBk{k}*vi));    % ignore small imaginary part
                    end
                    gainStatDesi = obj.Ph(n)*gainStatDesi;
                    obj.gainStatDes(n) = gainStatDesi;
                    
                    gainStatInti = 0;
                    QAmk = obj.oneSideCor(obj.thetA{m}, obj.del{m}, obj.nrx);
                    %QBmk = obj.oneSideCor(obj.thetB{m}, obj.del{m}, obj.ntx);
                    for k = 1:obj.ncluster(m)
                        gainStatInti = gainStatInti + real(obj.phi{m}(k)*(ui'*(QAmk{k})*ui)*txGain);
                        %gainStatInti = gainStatInti + real(obj.phi{m}(k)*(ui'*(QAmk{k})*ui)*obj.ntx);
                        %gainStatInti = gainStatInti + real(obj.phi{m}(k)*(ui'*(QAmk{k})*ui)*(vi'*QBmk{k}*vi));    % ignore small imaginary part
                    end
                    gainStatInti = obj.Ph(m)*gainStatInti;
                    obj.gainStatInt(n) = gainStatInti;
                end
                
                
                
                
                
                
                % second model computations
                
%                 if ~obj.skipSecComp
%                     QAkcon = QAk;
%                     for k = 1 : length(QAk)
%                         QAkcon{k} = conj(QAk{k});
%                     end
%                     
%                     
%                     SAs = obj.covMats(QBk, obj.thetB{n}, obj.thetA{n}, obj.del{n}, obj.ntx, obj.nrx);
%                     SBs = obj.covMats(QAk, obj.thetA{n}, obj.thetB{n}, obj.del{n}, obj.nrx, obj.ntx);
%                     RAQB = zeros(obj.nrx);
%                     RBQA = zeros(obj.ntx);
%                     phiKM = obj.phi{n}*(obj.phi{n})';
%                     for k = 1 : obj.ncluster(n)
%                         for m = 1 : obj.ncluster(n)
%                             RAQB = RAQB + SAs{k,m}*phiKM(k,m);
%                             RBQA = RBQA + SBs{k,m}*phiKM(k,m);
%                         end
%                     end
%                     [UA,~] = eig(RAQB);
%                     [UB,~] = eig(RBQA);
%                     
%                     W = zeros(obj.ntx,obj.nrx);
%                     for r = 1 : obj.nrx
%                         for t = 1 : obj.ntx
%                             for k = 1 : obj.ncluster(n)
%                                 W(t,r) = W(t,r) + obj.phi{n}(k)*real(((UA(:,r))'*(QAk{k})*UA(:,r))*((UB(:,t))'*QBk{k}*UB(:,t)));
%                             end
%                         end
%                     end
%                     
%                     [tempMax, tempI] = max(W);
%                     [gaini, rhat] = max(tempMax);
%                     that = tempI(rhat);
%                     ui = UA(:,rhat);
%                     vi = UB(:,that);
%                     obj.gainStatDes2(n) = gaini;
%                     obj.u2{n} = ui;
%                     obj.v2{n} = vi;
%                 end
%                 
                %---------------compute coherent gain begin ----------------
                
                
%                 % Generate random angles for each lobe
%                 angrx = (rand(nsub,obj.ncluster(n))-0.5)*diag(2*obj.del{n}) + repmat(obj.thetA{n}',nsub,1);
%                 %angtx = (rand(nsub,obj.ncluster(n))-0.5)*diag(2*obj.del{n}) + repmat(obj.thetB{n}',nsub,1);
%                 
%                 angrx2 = (rand(nsub,obj.ncluster(m))-0.5)*diag(2*obj.del{m}) + repmat(obj.thetA{m}',nsub,1);
%                 %angtx2 = (rand(nsub,obj.ncluster(m))-0.5)*diag(2*obj.del{m}) + repmat(obj.thetB{m}',nsub,1);
%                 
%                 % Vectorize angles and powers of all the subpaths
%                 p = repmat(obj.phi{n}'/nsub,nsub,1);
%                 p = p(:);
%                 angrx = angrx(:);
%                 %angtx = angtx(:);
%                 npath = nsub*obj.ncluster(n);
%                 
%                 p2 = repmat(obj.phi{m}'/nsub,nsub,1);
%                 p2 = p2(:);
%                 angrx2 = angrx2(:);
%                 %angtx2 = angtx2(:);
%                 npath2 = nsub*obj.ncluster(m);
%                 
%                 % Compute angular response
%                 erx = exp(2*pi*1i*cos(angrx)*(0:obj.nrx-1)*0.5);
%                 %etx = exp(2*pi*1i*cos(angtx)*(0:obj.ntx-1)*0.5);
%                 etx = ones(size(erx,1),obj.ntx);
%                 
%                 
%                 erx2 = exp(2*pi*1i*cos(angrx2)*(0:obj.nrx-1)*0.5);
%                 %etx2 = exp(2*pi*1i*cos(angtx2)*(0:obj.ntx-1)*0.5);
%                 etx2 = ones(size(erx2,1),obj.ntx);
%                 
%                 % Generate random fading realizations while keeping the angles fixed
%                 % For each realization, compute the optimal BF gain
%                 nt = 500;
%                 if ~obj.skipSecComp
%                     gainDes2t = zeros(nt,1);
%                 end
%                 if ~obj.minComp
%                     gainCohDest = zeros(nt,1);
%                 end
%                 gainIntt = zeros(nt,1);
%                 gainDest = zeros(nt,1);
%                 for t = 1:nt
%                     g = (randn(npath,1)+1i*randn(npath,1)).*sqrt(p/2);
%                     g = repmat(g,1,obj.nrx);
%                     H = (g.*erx)'*etx;
%                     [U,S,V] = svd(H,'econ');
%                     [s,ind]=max(diag(S));
%                     if ~obj.minComp
%                         gainCohDest(t) = s^2;
%                     end
%                     ut = U(:,ind);
%                     %vt = V(:,ind);
%                     tempv = (1/sqrt(obj.ntx))*ones(obj.ntx,1);
%                     
%                     %gainDest(t) = abs((obj.u{n})'*H*obj.v{n}).^2;
%                     gainDest(t) = abs((obj.u{n})'*H*tempv).^2;
%                     
%                     if ~obj.skipSecComp
%                         %gainDes2t(t) = abs((obj.u2{n})'*H*obj.v2{n}).^2;
%                         gainDes2t(t) = abs((obj.u2{n})'*H*tempv).^2;
%                     end
%                     
%                     g = (randn(npath2,1)+1i*randn(npath2,1)).*sqrt(p2/2);
%                     g = repmat(g,1,obj.nrx);
%                     H = (g.*erx2)'*etx2;
%                     gainIntt(t) = abs((ut'*H*tempv)).^2;
%                     %gainIntt(t) = abs((ut'*H*vt)).^2;
%                 end
%                 if ~obj.minComp
%                     obj.gainCohDes(n) = mean(gainCohDest);
%                 end
%                 obj.gainInt(n) = mean(gainIntt);
%                 obj.gainDes(n) = mean(gainDest);
%                 if ~obj.skipSecComp
%                     obj.gainDes2(n) = mean(gainDes2t);
%                 end
                
                %---------------compute coherent gain end ----------------
                
%                 if ~obj.minComp
%                     disp([num2str(n) '/' num2str(obj.nsample)]);
%                 end
            end
        end
        
        function plotData(obj)
            normalization = false;
            if normalization
                divdr = (obj.ntx*obj.nrx);
            else
                divdr = 1;
            end
            
            figure;
            tot = obj.nsample;
            plot(10*log10(sort(obj.gainDes)/divdr), (1:tot)/tot, 'Linewidth', 2, 'Color', [0 0 1]);
            hold on;
            plot(10*log10(sort(obj.gainInt)/divdr), (1:tot)/tot, 'Linewidth', 2, 'Color', [1 0 0]);
            legStr = {'ensemble average gain with kronecker model', 'ensemble average interference with kronecker model'};
            if ~obj.minComp
                plot(10*log10(sort(obj.eigProdKron)/divdr), (1:tot)/tot, 'Linewidth', 2, 'Color', [0 0.5 0.5]);
                plot(10*log10(sort(obj.gainStatDes)/divdr), (1:tot)/tot, 'Linewidth', 2, 'Color', [0 1 1]);
                plot(10*log10(sort(obj.gainCohDes)/divdr), (1:tot)/tot, 'Linewidth', 2, 'Color', [0 0 0]);
                plot(10*log10(sort(obj.gainStatInt)/divdr), (1:tot)/tot, 'Linewidth', 2, 'Color', [1 0 1]);
                legStr = [legStr, {'kronecker eigenvalue product', 'stat avg gain with kronecker model', 'coherent gain', 'stat avg interference with kron'}];
            end
            if ~obj.skipSecComp
                plot(10*log10(sort(obj.gainStatDes2)/divdr), (1:tot)/tot, 'Linewidth', 2, 'Color', [0 1 0.5]);
                plot(10*log10(sort(obj.gainDes2)/divdr), (1:tot)/tot, 'Linewidth', 2, 'Color', [0 1 0]);
                legStr = [legStr, {'stat avg gain with new model', 'ensemble avg with new model'}];
            end
            grid on;
            title(['Gain geometry for ' num2str(obj.ntx) ' Tx and ' num2str(obj.nrx) ' Rx antennas']);
            xlabel('Gain (dB)');
            legend(legStr, 'Location', 'best');
            
        end
        
        function [Qs] = oneSideCor(obj,thet, del, antsz)
            K = length(thet);
            Qs = cell(K,1);
            for k = 1 : K
                Qs{k}=obj.intCor(thet(k),del(k),antsz);
            end
        end
        
        function [Qk] = intCor(obj,thet,del,antsz)
            
            delAnt = 0.5;  % antenna spacing / lambda
            
            nstep = 1000;
            
            stepsize = 2*del/nstep;
            x = cos(thet-del:stepsize:thet+del);
            pow = (0:antsz-1)';
            pow = -1i*2*pi*delAnt*pow;
            
            e1 = exp(pow*x(1));
            
            eN = exp(pow*x(end));
            
            fx = toeplitz(e1') + toeplitz(eN');
            for k = 1:nstep-1
                e = 2*exp(pow*x(k+1));
                fx = fx + toeplitz(e');
            end
            fx = fx/(2*nstep);
            
            Qk = fx;
        end
        
        function [Ss] = covMats(obj, Q, thet, thet2, del, antsz, antsz2)
            K = length(thet);
            Ss = cell(K,K);
            for k = 1 : K
                for m = 1 : K
                    Ss{k,m}=obj.intCov(Q{k}, thet(m), thet2(m), del(m), antsz, antsz2);
                end
            end
            
        end
        
        function [Smk] = intCov(obj, Q, thet, thet2, del, antsz, antsz2)
            delAnt = 0.5;  % antenna spacing / lambda
            
            nstep = 1000;
            
            stepsize = 2*del/nstep;
            x = cos(thet-del:stepsize:thet+del);
            y = cos(thet2-del:stepsize:thet2+del);
            powx = -1i*2*pi*delAnt*(0:antsz-1)';
            powy = -1i*2*pi*delAnt*(0:antsz2-1)';
            
            ex1 = exp(powx*x(1));
            ey1 = exp(powy*y(1));
            
            exN = exp(powx*x(end));
            eyN = exp(powy*y(end));
            
            fx = ey1*(ex1'*Q*ex1)*ey1'+eyN*(exN'*Q*exN)*eyN';
            for k = 1:nstep-1
                exk = exp(powx*x(k+1));
                eyk = exp(powy*y(k+1));
                fx = fx + eyk*(2*(exk'*Q*exk))*eyk';
            end
            fx = fx/(2*nstep);
            
            Smk = fx;
        end

%         function J = getImage(obj, phi, thet, del, u)
%             antsize = length(u);
%             lam_real = 300/obj.fc;   % mm
%             antSep_real = lam_real/2; % mm
%             ant_sep = 3; %px
%             pixlen = antSep_real/ant_sep; %mm per pixel
%             
%             lam = 2*ant_sep; %px
%             N = 777;
%             % att=2e-4;
%             att = 4.58;
%             
%             center = round((N-1)/2);
%             ypos = center;
%             %phi = pi/5;
%             pos = zeros(antsize,2);
%             pos(:,2)=ypos;
%             offset = round((N-((antsize-1)*(ant_sep+1)+1))/2);
%             for k = 1:antsize
%                 pos(k,1) = offset+(k-1)*(ant_sep+1)+1;
%             end
%             
%             %phase = phi*[1:nrx];
%             
%             
%             I = zeros(N);
%             for k = 1 : antsize
%                 for x = 1 : N
%                     for y = 1 : N
%                         dist = (x-pos(k,1))+1i*(y-pos(k,2));
%                         absd = abs(dist);
%                         I(y,x)=I(y,x)+absd^-att*exp(1i*(2*pi*absd/lam))*u(k);
%                     end
%                 end
%             end
%             I = abs(I);
%             Imin = nanmin(I(:));
%             for k = 1:antsize
%                 I(pos(k,2),pos(k,1))= Imin;
%             end
%             Imax = nanmax(I(:));
%             I = (I - Imin)/(Imax-Imin);
%             
%             J = zeros(N,N,3);
%             J(:,:,1) = I;
%             J(:,:,2) = I;
%             J(:,:,3) = I;
%             
%             for k = 1:antsize
%                 J(pos(k,2)+[-1 0 1],pos(k,1)+[-1 0 1],1)= 1;
%                 J(pos(k,2)+[-1 0 1],pos(k,1)+[-1 0 1],2)= 0;
%                 J(pos(k,2)+[-1 0 1],pos(k,1)+[-1 0 1],3)= 0;
%             end
%             
%             cp = [center,center];
%             %I2 = zeros(N);
%             minrad = antsize*ant_sep;
%             maxrad = round(N/3);
%             for x = 1 : N
%                 for y = 1 : N
%                     dist = (x-cp(1))-1i*(y-cp(2));
%                     absd = abs(dist);
%                     ang = angle(dist);
%                     if ang<0
%                         ang = ang+2*pi;
%                     end
%                     for k = 1 : length(phi)
%                         tpd = thet(k)+del(k);
%                         tmd = thet(k)-del(k);
%                         if ((ang<tpd)&&(ang>tmd))||...
%                                 ((ang+2*pi<tpd)&&(ang+2*pi>tmd))||...
%                                 ((ang-2*pi<tpd)&&(ang-2*pi>tmd))
%                             if(absd<phi(k)*maxrad+1)&&(absd>phi(k)*maxrad-1)
%                                 J(y,x,:) = [0 1 0];
%                             end
%                         end
%                         if (absd<phi(k)*maxrad)&&(absd>minrad)
%                             if ((ang<tpd+0.02)&&(ang>tpd-0.02))||((ang<tmd+0.02)&&(ang>tmd-0.02))||...
%                                     ((ang+2*pi<tpd+0.02)&&(ang+2*pi>tpd-0.02))||((ang+2*pi<tmd+0.02)&&(ang+2*pi>tmd-0.02))||...
%                                     ((ang-2*pi<tpd+0.02)&&(ang-2*pi>tpd-0.02))||((ang-2*pi<tmd+0.02)&&(ang-2*pi>tmd-0.02))
%                                 J(y,x,:) = [0 1 0];
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%         
%         function seeGainRms(obj)
%             diff = obj.eigProdKron - real(obj.gainStatDes);
%             maxK = max(obj.ncluster);
%             gainrms = zeros(1,maxK);
%             errMax = zeros(1,maxK);
%             clusterOcc = zeros(1,maxK);
%             for k = 1 : maxK
%                 temp = diff(find(obj.ncluster==k));
%                 if ~isempty(temp)
%                     clusterOcc(k) = length(temp);
%                     gainrms(k) = sqrt(sum(temp.^2)/clusterOcc(k));
%                     errMax(1,k) = max(abs(temp));
%                 end
%             end
%             figure;
%             [AX,H1,H2] = plotyy(1:maxK,errMax,1:maxK,clusterOcc,'bar','bar');
%             set(get(AX(1),'Ylabel'),'String','rms error');
%             set(get(AX(2),'Ylabel'),'String','occurence');
%             set(H1, 'BarWidth', 0.6, 'FaceColor', [0 0 1], 'DisplayName', 'Maximum Absolute Error');
%             set(H2, 'BarWidth', 0.1, 'FaceColor', [1 0 0], 'DisplayName', 'Number of Occurrence of k cluster');
%             hold on;
%             bar(gainrms, 'BarWidth', 0.3, 'FaceColor', [0 1 1], 'DisplayName', 'RMS Error between Kronecker and Actual Gain');
%             legend('show');
%             xlabel('# clusters');
%         end
    end
end