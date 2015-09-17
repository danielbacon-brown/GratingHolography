%Runs the built-in Matlab GA

function  RunGA_addCircles( GAoptions )

[chromosome, fitness, reason, output, pop, scores] = ga(GAoptions.GAproblem);
save([GAoptions.dir,GAoptions.GARecordFileBase,int2str(GAoptions.currentGAnum)],'chromosome','fitness','reason','output','pop','scores');

end