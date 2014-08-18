%%Driver

n = 16;
para = ParaClass([n,n,n]);
% surveyDesign = struct('type','monopole');
surveyDesign = struct('type','perms','padding',[1,1],'skip',3,'plotIt',true);
para = para.createSurvey(surveyDesign);
figure
modelDesign = struct('type','eh','position','middle','size',[6 6 10],'values',log([1E-2,1E-1]),'plotIt',true);
para = para.createModel(modelDesign);
para.mref = para.m.*0+modelDesign.values(1);
dataDesign = struct('type','normal','error',0,'outliers',0);
para = para.createData(dataDesign,0.02,1E-2);

para.beta = 10^1;
para.GRADw = para.createGRADw(1,1,1,1,0);


% para.objFnc = 'l2';para.beta = 1e4;para.GRADw = para.createGRADw(1,1,1,2,pi);
% [x, info] = minimize( para, para.mref, 1E-3, 10,struct('plotIt',true));