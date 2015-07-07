
% script to discover MMCU from dataset 1

%% load vertices
subList_HVU = importdata ('subjects_harvard.list');
subList_SWU = importdata ('subjects_swu.list');

subList = [subList_HVU; subList_SWU]; 
siteIdx = [ones(length(subList_HVU), 1); 3*ones(length(subList_SWU),1)];

metrics = {'alff', 'dc', 'ec', 'reho', 'reho2', 'falff',  'lgi',  'sulc', 'thickness','meancurv', 'area',  'volume'};

sites = {'HVU', 'SWU'};
nMetrics = length (metrics);


% load mixing time series
st = 1;
stNm = sites{st};
list = ['subList_', stNm],
site_code = st,

metricMix = [];
obj.setup.subNum = nMetrics;
obj.setup.trial  = ones(nMetrics, 1);
obj.result.trialTab = zeros(nMetrics,3); 

for m = 1:nMetrics
	metricNm = metrics{m},
	fn = sprintf ('../singleMetricICA_%s/%s_%s.ica/melodic_mix', stNm, stNm, metricNm);
	tmp = load (fn);
	metricMix = [metricMix, tmp];
	size(obj.result.trialTab),
	obj.result.trialTab(m,:) = [1, m, size(tmp, 2)];
end

% compute similarity matrix
obj.result.MICM = abs (corrcoef (metricMix));
for m = 1:nMetrics
	range = sum(obj.result.trialTab(1:m-1,3))+1:sum(obj.result.trialTab(1:m,3));
	obj.result.MICM(range, range) = 0;
end
MICM_orig = obj.result.MICM;

%% align components
obj = findComp_PR(obj, 0);

% examine repro
flag_genNull = 1;
numIter = 1000;


obj.result.MICM = MICM_orig;
obj = normByRow (obj);
obj.result.MICM = obj.result.MICM'+obj.result.MICM;
min_MICM = min (obj.result.MICM(:));
obj.result.MICM = obj.result.MICM - min_MICM;

if flag_genNull == 1
    null_interSb_reproMap = generateNullRepro (obj, numIter);
end

% [obj.result.beta_rank_subjLoad, obj.result.sig_subjLoad, obj.result.subjLoad,obj.result.null_subjLoad] = examSigSubjLoad (obj, null_interSb_reproMap);
obj = examConnections (obj, null_interSb_reproMap);


%% compute similarity between matched metrics using correlation
base_pos = cumsum(obj.result.trialTab(:,3));
base_pos = [0; base_pos(1:end-1) ];
numCp = size(obj.result.foundComp, 1);
cp_cc_mtx = zeros (nMetrics, nMetrics, numCp);
for cp = 1:numCp
    cps = obj.result.foundComp(cp, :);
    idx_rm = find (cps == 0);
    cps(idx_rm) = 1;
    idx = base_pos' + cps;
    cp_cc_mtx(:,:,cp) = MICM_orig(idx, idx);
    cp_cc_mtx(idx_rm, :, cp) = 0;
    cp_cc_mtx(:, idx_rm, cp) = 0;
end
    
