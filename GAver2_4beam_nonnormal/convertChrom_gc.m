function varargout = convertChrom_gc(chromosome, sectionLengths)
%Converts the chromosome or chromosome section into array of the corresponding fractional values 0<=x<1
%Uses graycode

u=0;
cellout = cell(size(sectionLengths,1)); %preallocate output 

for iter = 1:max(size(sectionLengths)) %number of sections   
    cellout{iter} = (gc2dec( chromosome(u+1:u+sectionLengths(iter)) )+1)/(2^sectionLengths(iter)); %Selects section and converts to fraction %not allowed to be 0
    u = u + sectionLengths(iter);  %Moves pointer to next section
end

%if only one argument is asked for, will return as a vector
if nargout == 1
    varargout = {cell2mat(cellout)};
else
    varargout = cellout;
end


end