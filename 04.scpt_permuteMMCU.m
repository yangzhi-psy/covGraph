
%% scpt for permuting SCMs to generate fake MMCUs and compute their ICCs across sites

%% load spatial maps from HVU (templateMap)
metrics = {'alff', 'dc', 'ec', 'reho', 'reho2', 'falff',  'lgi',  'sulc', 'thickness','meancurv', 'area',  'volume'};

 stNm = 'HVU'
 
 templateMap = cell(length(metrics), 1);
 for m = 1:length(metrics)
     metricNm = metrics{m},
 	fn = sprintf ('../singleMetricICA_%s/%s_%s.ica/melodic_IC.nii.gz', stNm, stNm, metricNm);
     fn = gunzip (fn);
 	nii = load_nii (cell2mat(fn));
     
     sz = size (nii.img);
     map2D = zeros(sz(1)*sz(2), sz(3));
     for c = 1:sz(3)
         map2D(:, c) = reshape (nii.img(:,:,c), [sz(1)*sz(2), 1]);
     end
     
     templateMap(m) = {map2D};
 end
 
 %% spatial regression to get mixing time series in SWU
 load cp_cc_mtx_HVU.mat
 nCp = length (obj.result.foundRepro);
 
 % load single subject data
 load zeroVertIdx
 subList = importdata ('subjects_swu.list');
 siteIdx = 3*ones(length(subList),1);
 
 nIt = 5000
 
 mixT = zeros (length(subList), length(metrics), nIt);
 for m = 1:length(metrics)
     m,
     
     vertexData = fnc_collectVertices (subList, siteIdx, metricNm);
     vertexData = vertexData';
     vertexData(zeroVertIdx, :) = [];

     rndIdx(:, m) = randi([0,max(obj.result.foundComp(:,m))], nIt,1);
     
     for cp = 1:nIt
         idx_IC = rndIdx(cp, m);
         if idx_IC == 0
             mixT(:, m, cp) = zeros (size(vertexData, 2), 1);
         else
             t_map = templateMap{m}(:, idx_IC);
             mixT(:,m, cp) = corr(t_map, vertexData);
         end
     end
 end
 save mixT_SWU_permuteMMCU_.mat mixT rndIdx

%% compute correlation matrices for faked components

load mix_container.mat % mix series for HVU

load mixT_SWU_permuteMMCU_.mat

metrics = {'alff', 'dc', 'ec', 'reho', 'reho2', 'falff',  'lgi',  'sulc', 'thickness','meancurv', 'area',  'volume'};
nIt = 5000
% fake corr mat for HVU
mixT_HVU = zeros(184,12,nIt);
for m = 1:length(metrics)
    m,
    tcs = mix_container{m};
    for cp = 1:nIt
        idx_IC = rndIdx(cp, m);
        if idx_IC ~= 0
            mixT_HVU(:, m, cp) = tcs(:, idx_IC);
        end
    end
end


fake_cc_mtx_SWU = zeros(12, 12, nIt);
fake_cc_mtx_HVU = zeros(12, 12, nIt);

for cp = 1:nIt
    fake_cc_mtx_SWU(:,:,cp) = abs(corrcoef(mixT(:,:,cp))) - eye(12);
    fake_cc_mtx_HVU(:,:,cp) = abs(corrcoef(mixT_HVU(:,:,cp))) - eye(12);
end

save fake_cc_mtx.mat fake_cc_mtx_SWU fake_cc_mtx_HVU


%% compute ICCs
fake_cc_mtx_SWU(isnan(fake_cc_mtx_SWU)) = 0;
fake_cc_mtx_HVU(isnan(fake_cc_mtx_HVU)) = 0;


null_crsST_ICC = zeros (nIt, 1);
for cp = 1:nIt
    scan1 = squareform(fake_cc_mtx_HVU(:,:,cp), 'tovector');
    scan2 = squareform(fake_cc_mtx_SWU(:,:,cp), 'tovector');
    null_crsST_ICC(cp) = ICC([scan1', scan2'], '1', '1');
end

save null_crsST_ICC_permuteMMCU.mat null_crsST_ICC
