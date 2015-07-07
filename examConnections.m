
function obj = examConnections (obj, null_interSb_reproMap)

sz = size (null_interSb_reproMap);

null_conn = zeros (obj.setup.subNum, obj.setup.subNum, sz(2));
for i = 1:sz(2)
    tmp_map = squareform (null_interSb_reproMap(:,i), 'tomatrix');
    null_conn(:,:,i) = tmp_map;
end

numComp = length (obj.result.foundRepro);
p_conn = zeros(obj.setup.subNum, obj.setup.subNum, numComp);

for i = 1:numComp
    
    fprintf ('computing component %d of %d\n', i, numComp);
    reproMap = cell2mat(obj.result.foundRepro(i));
    reproMap = full (reproMap+reproMap');
    similarity = sepIntra_Inter (obj, [], reproMap);
    p_conn(:,:,i) = sum(repmat (similarity, [1,1,sz(2)]) < null_conn, 3)./sz(2);
end

obj.result.p_conn = p_conn;


    
    