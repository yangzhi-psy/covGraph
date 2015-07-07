function vertexData = fnc_collectVertices (subList, siteIdx, metricNm)
nVertices_lh = 10242;
nVertices_rh = 10242;


vertexData = zeros (length(subList), nVertices_lh+nVertices_rh);

if strcmp (metricNm, 'reho')
	sm = 0;
elseif strcmp (metricNm, 'reho2')
	sm = 0;
else
	sm = 4;
end

for i = 1:length (subList)
	sbNm = subList{i};
	switch siteIdx(i)
		case 1
			siteNm = 'harvard_buckner';
		case 2
			siteNm = 'nki_lifespan';
			sbNm = sbNm(4:end);
		case 3
			siteNm = 'swu_qiu';
	end

	fn = sprintf ('../%s/%s/fsaverage5/lh.%s.sm%d.nii.gz', siteNm, metricNm, sbNm,sm);
	tmphead = load_nifti(fn);
	vertexData(i, 1:nVertices_lh) = tmphead.vol';

	fn = sprintf ('../%s/%s/fsaverage5/rh.%s.sm%d.nii.gz', siteNm, metricNm, sbNm,sm);
	tmphead = load_nifti(fn);
	vertexData(i, nVertices_lh+1:nVertices_lh+nVertices_rh) = tmphead.vol';
end
