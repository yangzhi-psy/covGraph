% script to validate MMCUs in dataset 2.
% the MMCUs found in dataset1 are used as spatial templates
% spatial regression is used to compute mixing course in dataset2
% similarity matrices between mixing courses obtained in dataset2 are then computed
% ICC are computed between the similarity matrics obtained in dataset1 and dataset2.

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

%% spatial regression
load cp_cc_mtx_HVU.mat
nCp = length (obj.result.foundRepro);

% load single subject data
load zeroVertIdx
subList = importdata ('subjects_swu.list');
siteIdx = 3*ones(length(subList),1);


mixT = zeros (length(subList), length(metrics), nCp);
for m = 1:length(metrics)
    m,
    
    vertexData = fnc_collectVertices (subList, siteIdx, metricNm);
    vertexData = vertexData';
    vertexData(zeroVertIdx, :) = [];
    
    for cp = 1:nCp
        idx_IC = obj.result.foundComp(cp, m);
        if idx_IC == 0
            mixT(:, m, cp) = zeros (size(vertexData, 2), 1);
        else
            t_map = templateMap{m}(:, idx_IC);
            mixT(:,m, cp) = corr(t_map, vertexData);
        end
    end
end
save mixT_SWU.mat mixT

%% compare mixPatterns between template and regression datasets
load mixT_SWU.mat
load cp_cc_mtx_HVU.mat
metrics = {'alff', 'dc', 'ec', 'reho', 'reho2', 'falff',  'lgi',  'sulc', 'thickness','meancurv', 'area',  'volume'};

nCp = length (obj.result.foundRepro);

cp_cc_mtx_reg = zeros (size(cp_cc_mtx));
for cp = 1:nCp
    cp_cc_mtx_reg(:,:, cp) = abs(corrcoef (mixT(:,:,cp))) - eye(length(metrics));
end
idx_rm = isnan(cp_cc_mtx_reg);
cp_cc_mtx_reg(isnan(cp_cc_mtx_reg)) = 0;
%  figure, for i = 1:40,imagesc([cp_cc_mtx(:,:,i), cp_cc_mtx_reg(:,:,i)]);i,pause;end    

% convert cp_cc_mtx to disPat
disPat_HVU = fnc_ccMtx2disPat (cp_cc_mtx);
disPat_SWU = fnc_ccMtx2disPat (cp_cc_mtx_reg);

% cross-trial correlation
crsST_cc_mtx = corrcoef ([disPat_HVU, disPat_SWU]);

for cp = 1:nCp
    repro(cp) = crsST_cc_mtx(cp, cp+size(cp_cc_mtx, 3));
    df(cp) = nnz(disPat_SWU(:,cp)-1);
    [t_repro(cp), p_repro(cp)] = r2t(repro(cp), df(cp));
end

t_repro(df<2) = 0;
p_repro(df<2) = 1;

% ICC and Kendall's W
for cp = 1:nCp
	if t_repro(cp) ~= 0
%	tmp = std (disPat_HVU(:, cp) - disPat_SWU(:, cp));
%	typicalErr(cp) = tmp/sqrt(2);
    tmp_idx = find (disPat_HVU(:, cp) ~= 1);
	crsST_ICC(cp) = ICC([disPat_HVU(tmp_idx, cp), disPat_SWU(tmp_idx, cp)], '1', '1');
    crsST_W(cp) = IPN_kendallW([disPat_HVU(tmp_idx, cp), disPat_SWU(tmp_idx, cp)], 0);
	end
end

save crsST_cc_HVU_SWU.mat repro crsST_cc_mtx cp_cc_mtx cp_cc_mtx_reg t_repro p_repro crsST_ICC;
