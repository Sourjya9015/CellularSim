

AntenNum = {[4,4], [4,8], [8, 4], [8,8], [16,16]};
%AntenNum = {[4,4]};

for i = 1:length(AntenNum)
    ntx = AntenNum{i}(1);
    nrx = AntenNum{i}(2);
    b = BeamFormingGain(ntx,nrx);
    gainMat = zeros(length(b.gainStatDes),2);
    %gainMat(:,1) = 10*log10(b.gainStatDes);
    gainMat(:,1) = b.gainStatDes;
    %gainMat(:,2) = 10*log10(b.gainStatInt);
    gainMat(:,2) = b.gainStatInt;
    gainMat = sort(gainMat);
    str = sprintf('bfGain%dX%d',ntx,nrx);
    cmd = ['save data/', str, ' gainMat'];
    eval(cmd);
end
