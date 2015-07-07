% scpt to compute overall co-variance between 12 metrics

load crsST_cc_HVU_SWU.mat
load cp_cc_mtx_HVU.mat

toshow_idx = find (crsST_ICC>0.6);
[sorted_crsST_ICC, sort_idx] = sort (crsST_ICC(toshow_idx), 'descend');
sorted_crsST_ICC,
toshow_idx = toshow_idx(sort_idx);


%%
toaverage = cp_cc_mtx (:, :, toshow_idx(1:end)); %1:99); %

r_mean_highRep =zeros(12);
for i = 1:12
    for j = i+1:12
        r_mean_highRep(i, j) = mean(r2z(nonzeros(toaverage(i, j,:))));
    end
end

r_mean_highRep = r_mean_highRep+r_mean_highRep';
Z = linkage(squareform(2-r_mean_highRep-2*eye(12), 'tovector'));

% figure, dendrogram(Z, 0)
figure, imagesc(r_mean_highRep)

figure, [h, t, outperm] = dendrogram(Z, 0);
a = r_mean_highRep;
a = a(outperm, :);
a = a(:, outperm);
figure, imagesc(a);
