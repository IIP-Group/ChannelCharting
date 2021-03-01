% =========================================================================
% -- Functions to calculate trustworthiness (TW) and continuity (CT)
% -------------------------------------------------------------------------
% -- (c) 2016-2021 Christoph Studer, Emre Gonultas, and Said Medjkouh
% -- e-mail: studer@ethz.ch, eg566@cornell.edu, sm2685@cornell.edu
% =========================================================================

function TW_CT(location,mappedX,par)

index = 1:par.U;
knMax = round(0.05*par.U);
kn= ceil(linspace(1,knMax,6));

N = length(index); %total number of users

disp('Calculating pairwise distances in both domains');

d = zeros(size(location,2),size(location,2));
for qq=1:size(location,2)
    d(:,qq)= sqrt(  sum( bsxfun(@minus,location,location(:,qq)).^2)  );
end

d_cov=zeros(size(mappedX,1),size(mappedX,1));
for qq=1:size(mappedX,1)
    d_cov(:,qq)= sqrt(  sum( bsxfun(@minus,mappedX',mappedX(qq,:)').^2)  );
end

d(d==0)=inf;
d_cov(d_cov==0)=inf;

[~,inds]=sort(d); %inds(:,i) contains i-neighbors (real space) indexes in order
[~,inds_cov]=sort(d_cov); %inds_cov(:,i) contains i-neighbors (covMat space) indexes in order

tw =zeros(1,length(kn)); %trustworthiness
ct =zeros(1,length(kn)); %continuity

idx = 0;
disp('Calculating TW and CT...');
for k=kn
    
    idx = idx +1;
    
    fprintf('Calculating for %d neighbors...\n',k);
    
    for i=1:length(location) %for all points in x
        % Uki: indices of elements that are among k-NN of i (in covMat) but not in real locations
        Uki = intersect(inds_cov(1:k,i),inds(k+1:end,i));
        
        % find the ranks of Uki elements in the real locations rank
        [~,rankUki] = ismember(Uki,inds(:,i));
        twi = sum(rankUki - k);
        tw(idx) = tw(idx) + twi;
        
        % Vki: indices of elements that are among k-NN of i (in covMat) but not in real locations
        Vki = intersect(inds_cov(k+1:end,i),inds(1:k,i));
        
        %find the ranks of Uki elements in the real locations rank
        [~,rankVki] = ismember(Vki,inds_cov(:,i));
        cti = sum(rankVki - k);
        ct(idx) = ct(idx) + cti;
        
    end
end

% compute TW and CT
TW = 1 - 2*tw./(N*kn.*(2*N-3*kn-1));
CT = 1 - 2*ct./(N*kn.*(2*N-3*kn-1));

strFigTrust = sprintf('%s U%d B%d ',par.nameParam, par.U,par.B);

% plot results
h = figure();
plot(kn,TW,'bo-');
hold on
plot(kn,CT,'rs-');
hold off
xlabel('k-nearest neighbors')
ylabel('TW and CT');
title(strrep(strFigTrust,'_',' '))
legend('TW','CT','Location','southwest')
grid on
axis([0 max(kn) 0 1])

% save figure (in color and with a reasonable bounding box)
strFigTrust = sprintf('TW_CT_%s_U%d_B%d',par.nameParam, par.U,par.B);
print(['./output/',strFigTrust,'.png'],'-dpng')

end