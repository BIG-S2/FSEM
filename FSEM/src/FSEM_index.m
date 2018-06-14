function [MZtp1, MZtp2, DZtp1, DZtp2, MDti] = FSEM_index(famid, zygosity)
% FSEM_index is aimed to extract the corresponding 
%   MZ(DZ) twin pairs and twin individual's indices 

famid = sort(famid);

indt1 = ([famid(2:end);0] == famid);
indt2 = ([0;famid(1:(end-1))] == famid);

% they are all boolean vectors 
MZtp1 = (indt1 & zygosity == 0);
MZtp2 = (indt2 & zygosity == 0);
DZtp1 = (indt1 & zygosity == 1);
DZtp2 = (indt2 & zygosity == 1);
MDti = ~(MZtp1 | MZtp2 | DZtp1 | DZtp2);

end
