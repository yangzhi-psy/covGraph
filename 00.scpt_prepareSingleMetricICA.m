
% script to combine metric maps from all subjects for ICA.

%fs_home = '/opt/freesurfer'; fsaverage = 'fsaverage5';
%fannot = [fs_home '/subjects/' fsaverage '/label/lh.aparc.a2009s.annot'];
%vertices_lh = read_annotation(fannot);
%nVertices_lh = numel(vertices_lh);
%
%fannot = [fs_home '/subjects/' fsaverage '/label/rh.aparc.a2009s.annot'];
%vertices_rh = read_annotation(fannot);
%nVertices_rh = numel(vertices_rh);
%clear fannot

%% load vertices
subList_HVU = importdata ('subjects_harvard.list');
subList_SWU = importdata ('subjects_swu.list');

subList = [subList_HVU; subList_SWU]; 
siteIdx = [ones(length(subList_HVU), 1); 3*ones(length(subList_SWU),1)];

metrics = {'alff', 'area', 'dc', 'ec', 'falff', 'lgi', 'meancurv', 'reho', 'reho2', 'sulc', 'thickness', 'volume'};

allZeroVertIdx = [];
for m = [1:12] 
	metricNm = metrics{m},
	vertexData = fnc_collectVertices (subList, siteIdx, metricNm);
	% examine zero vertices that are consistent across subjects
	sumvert = sum(abs(vertexData), 1);
	allZeroVertIdx = [allZeroVertIdx, find(sumvert == 0)];
end

zeroVertIdx = unique(allZeroVertIdx);

% % load pre-stored zero vertices index
%load  zeroVertIdx.mat;
load nii_template.mat;

list = 'subList_HVU', % subList_SWU
site_code = 1, % 3

for m = 1:12
	metricNm = metrics{m},
	vertexData = fnc_collectVertices (eval(list), siteIdx(siteIdx==site_code), metricNm);
	tmp = vertexData';
	tmp(zeroVertIdx, :) = [];

	% scale
	gmin = min(tmp(:));
	gmax = max(tmp(:));
	tmp = (tmp - gmin)./(gmax - gmin)*10000;

	img = zeros (47,383,1, length(eval(list)));
	for i = 1:length(eval(list))
		img(:,:,1,i) = reshape (tmp(:,i), [47, 383]);
	end

	nii.img = img;
	nii.hdr.dime.dim = [4 47 383 1 length(eval(list)) 1 1 1];
	nii.hdr.glmax    = max(img(:));
	nii.hdr.glmin    = min(img(:));

	fn = sprintf ('../singleMetricICA_%s/%s_forICA_%s.nii', list(end-2:end), list(end-2:end), metricNm);
	save_nii (nii, fn);
end
