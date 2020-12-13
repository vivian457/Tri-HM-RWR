% nfold = 5: local 5fold CV 局部五折
% nfold = 1: novel disease prediction 全新疾病预测   
% clear;  fin = 'RWR_HM_DATA1-xj.mat';data=load(fin);
% clear;  fin = 'RWRHM_DATA1_Intersection-xj.mat';data=load(fin);
clear;  fin = 'RWR_HM_DATA1-update.mat';data=load(fin);
fdatestr =datestr(now,'yyyy.mmm.dd-HH.MM.SS');
prefix = 'local_5CV-DisByDis-'
rand('twister',mod(floor(now*8640000),2^31-1))
% % 
% % 原文默认参数 
% % options.alpha = 0.45;
% options.L1 = 4;
% options.r1 = 1;
% options.L2 = 1; 
% options.r2 = 0; 
options.alpha = 0.7; 
options.L1 = 1;   %对应 gene 相关参数
options.r1 = 2 ;   
options.L2 = 1 ;   %对应 miRNA 相关参数
options.r2 = 6; 
options.pro_jump = 0.5;  % 层间跳转概率   
options.NormalizationType = 'row'; 
% options.NormalizationType = 'col'; 
AsNewMiRNA = 0   ; 
knn=300       
% %  
Aprotein_disease     = sparse( data.hGDA   )   ;   % corresponding to miRNAs in this dataset     
BmiRNA_disease       =  ( data.hMDA   )  ;   [n_miRNA,n_disease ]= size( BmiRNA_disease); 
Cprotein_miRNA       = sparse( data.hMGA'  )  ; 
% % kp = sum(Cprotein_miRNA,2)+1; 
% % Aprotein_disease = Aprotein_disease(kp>0,:);
% % Cprotein_miRNA   = Cprotein_miRNA(kp>0,:);
% Cprotein_miRNA       = sparse( double( Aprotein_disease*BmiRNA_disease' >0 ) )  ; 
% % 
Protein_adjset.eye    = speye(  size(data.ggPathway)  );
% Protein_adjset.eye    = speye(  nnz(kp)  );
Protein_adjset.ggCoex_Interaction    =data.ggCoexpress;
Protein_adjset.ggPathway_Interaction = data.ggPathway; 
Protein_adjset.ggPPI_Interaction     =data.ggPPI; 
% % % Protein_adjset.ggSim_GDA = sparse( getAdjKnnColumns( getCosineSmilarityOfRowVectors( Aprotein_disease),  knn , 1, 1 ) ) ;         % GO_adjset corresponding to diseases in this dataset   
% % 
Disease_adjset.ddSS_394 = sparse( data.ddSS_394 );           % GO_adjset corresponding to diseases in this dataset   
% Disease_adjset.ddSSP = sparse( data.ddSSP );        % GO_adjset corresponding to diseases in this dataset   
% Disease_adjset.ddSemSIM_394 = sparse( getAdjKnnColumns( data.ddSemSIM_394,  knn , 1, 1 )  );    n_disease = size(Disease_adjset.ddSemSIM_394,1);       % GO_adjset corresponding to diseases in this dataset   
% Disease_adjset.ddSemSIM_GDA = getAdjKnnColumns( getCosineSmilarityOfRowVectors( Aprotein_disease'),  knn , 1, 1 )  ;         % GO_adjset corresponding to diseases in this dataset   
% Disease_adjset.ddSemSIM_MDA = getAdjKnnColumns( getCosineSmilarityOfRowVectors( BmiRNA_disease'),  knn , 1, 1 ) ;         % GO_adjset corresponding to diseases in this dataset   
% %  
miRNA_adjset.mmFS_678   = sparse( data.mmFS_678 )   ;  % Domain_adjset  corresponding to miRNAs in this dataset     
miRNA_adjset.mmSS_678   = sparse( data.mmSS_678 )   ;  % Domain_adjset  corresponding to miRNAs in this dataset     
% miRNA_adjset.mmFSP   = sparse( data.mmFSP )   ; % Domain_adjset  corresponding to miRNAs in this dataset     
% % miRNA_adjset.mmMirsim_678   = sparse( getAdjKnnColumns( data.mmMirsim_678,  knn , 1, 1 )  );     % Domain_adjset  corresponding to miRNAs in this dataset     
% miRNA_adjset.mmMirsim_GMA   = getAdjKnnColumns( getCosineSmilarityOfRowVectors(Cprotein_miRNA') ,  knn , 1, 1 )  ;     
% miRNA_adjset.mmMirsim_GMA2   = getAdjKnnColumns( getCosineSmilarityOfRowVectors(Cprotein_miRNA') ,  knn , 1, 1 )  ;     
% miRNA_adjset.mmMirsim_MDA    = getAdjKnnColumns( getCosineSmilarityOfRowVectors(BmiRNA_disease)  ,  knn , 1, 1 )  ;     
% % 
n_miRNAOfDiseases = sum(BmiRNA_disease,1);  %% number of related mirnas with each disease
ind_pos = find(BmiRNA_disease);
ind_neg = find(~BmiRNA_disease);
n_pos= length(ind_pos);
nfold    = 5;  
% nfold = n_pos
nCVTimes = 1; 
% %
% Indices = crossvalind('Kfold', n_pos, nfold); 
% EffDisIDset = find(n_miRNAOfDiseases>=nfold);  
EffDisIDset = find(n_miRNAOfDiseases>=21);
EffDisNames = data.diseaseList(EffDisIDset);
n_miRNAOfDiseases_EffDis=n_miRNAOfDiseases(EffDisIDset);  %% number of related mirnas with effective diseases
tbRes = table; 
for iDis=1:length(EffDisIDset)
    DisID = EffDisIDset(iDis);
    DisIDstr = ['DisInd:',num2str(DisID)]; 
    label = full( BmiRNA_disease(:,  DisID   ) ); 
    ind_pos = find(label);  n_pos= length(ind_pos);
    ind_neg = find(~label); n_neg= length(ind_neg);    
    AUCset = [] ;
    AUPRset = [] ;
	for iCV = 1:nCVTimes
        Indices = crossvalind('Kfold', n_pos, nfold);         
        Indices_neg = crossvalind('Kfold', n_neg, nfold);         
        for ifold=1:nfold
            ind_pos_test = ind_pos(Indices==ifold); 
            ind_neg_test = ind_neg(Indices_neg==ifold); 
            ind_test = [ind_neg_test; ind_pos_test ] ;  % 未知样本也被平均分为等份，抽出一份未知样本和一份正样本组成了一份测试样本  
% %           ind_test = [ind_neg; ind_pos_test ] ;    % 所有的未知样本作为负样本  
            label_train = label ; label_train(ind_test) = 0; 
            label_test  = label ; label_test(~ind_test) = 0 ; 
            BmiRNA_disease_train = BmiRNA_disease ; BmiRNA_disease_train(:,DisID) = label_train; 
% % % % %             BmiRNA_disease_test  = zeros( size(BmiRNA_disease) ) ; BmiRNA_disease_test(:,DisID) = label_test; 
            if  AsNewMiRNA
                BmiRNA_disease_train( ind_pos_test , :) = 0 ;    % 删除 测试miRNA的所有关联   
            end  
            % [kd,km] = gaussiansimilarity(BmiRNA_disease_train,n_disease,n_miRNA);  
            % Disease_adjset.ddSemSIM_MDA_pGIP = kd; %getAdjKnnColumns( kd , knn , 1, 1 ) ;         % GO_adjset corresponding to diseases in this dataset   
            % miRNA_adjset.mmMirsim_MDA_pGIP   =km; % getAdjKnnColumns( km  ,  knn , 1, 1 )  ;     
            Disease_adjset.ddSemSIM_MDA_cos = getAdjKnnColumns( getCosineSmilarityOfRowVectors( BmiRNA_disease_train'), knn , 1, 1 ) ;         % GO_adjset corresponding to diseases in this dataset   
            miRNA_adjset.mmMirsim_MDA_cos   = getAdjKnnColumns( getCosineSmilarityOfRowVectors(BmiRNA_disease_train)  ,  knn , 1, 1 )  ;     
            % %  [sd,sm] = integratedsimilarity(FS,FSP,SS,SSP,kd,km)   
            % [kd,km] = gaussiansimilarity(BmiRNA_disease_train,n_disease,n_miRNA);  
            %  [sd,sm] = integratedsimilarity(miRNA_adjset.mmFS_678,data.mmFSP,Disease_adjset.ddSS_394,data.ddSSP,kd,km);
            % % %  miRNA_adjset.mmFS_678_pGIP  = sm ;
            % % %  Disease_adjset.ddSS_394_pGIP = sd  ;
            %  miRNA_adjset.mmFS_678  = sm ;
            %  Disease_adjset.ddSS_394 = sd  ;
                [RmiRNA_dis, Rprotein_dis   , RDset, RPset, kp ]= A_ThrRW_Multiplex(Protein_adjset,Disease_adjset,miRNA_adjset,   Aprotein_disease,BmiRNA_disease_train,Cprotein_miRNA, options );
%                 [Rprotein_dis   ,RmiRNA_dis,  RDset, RPset, kp ]= A_ThrRW_Multiplex(miRNA_adjset,Disease_adjset,Protein_adjset,   BmiRNA_disease_train,Aprotein_disease,Cprotein_miRNA', options );
 
            [X,Y,T,AUC] = perfcurve(BmiRNA_disease(ind_test,DisID ),double(RmiRNA_dis(ind_test,DisID)),1) ;  
            [X_pr,Y_pr,tpr_pr,aupr_krr] = perfcurve( BmiRNA_disease(ind_test,DisID ),double(RmiRNA_dis(ind_test,DisID))  ,1, 'xCrit', 'reca', 'yCrit', 'prec');
            AUCset(ifold,iCV) = AUC ;
            AUPRset(ifold,iCV) = aupr_krr ;
        end
	end
    auc = mean(AUCset(:));
    std_auc=std(AUCset(:));
    auprc = mean(AUPRset(:));
    std_auprc=std(AUPRset(:)); 
%     [iDis n_miRNAOfDiseases_EffDis(iDis) auc std_auc  auprc  std_auprc ] 
%     tbRes.AUROC(iDis,1) = auc;
%     tbRes.AUPRC(iDis,1) = auprc;
    tbRes.AUROC(DisIDstr,1) = auc;
    tbRes.AUPRC(DisIDstr,1) = auprc;
    tbRes.std_auc(DisIDstr,1) = std_auc;
    tbRes.std_auprc(DisIDstr,1) = std_auprc;
    tbRes.relatedMirnaNum(DisIDstr,1) = n_miRNAOfDiseases_EffDis(iDis);
    tbRes
end
mean(tbRes{:,:},1) 

n_dis = size(tbRes,1); 
tbRes.DisName = EffDisNames ;
tbRes.AUROC('mean',1) = mean(tbRes.AUROC) ;  
tbRes.AUPRC('mean',1) = mean(tbRes.AUPRC);   
% tbRes
fout = [prefix,'-nfold=',num2str(nfold),'-nTimes=',num2str(nCVTimes),'-',fdatestr,'-','res-',fin,'.mat']
save(fout, 'tbRes',  'nfold',  'knn', 'options' );


