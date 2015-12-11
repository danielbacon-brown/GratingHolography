function [varargout] = splitChromosome(chromosome, sectionLengths)
%Splits the chromosome into seperate chromosomes based on the sectionLengths

u=0;
for iter = 1:max(size(sectionLengths)) %number of sections   
    varargout{iter} = chromosome( u+1: u+sectionLengths(iter) ); %Selects section and converts to fraction
    u = u + sectionLengths(iter);  %Moves pointer to next section
end


end