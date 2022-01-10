% Calculate functional connectivity (orthogonalized power envelope correlation) with Band-passed complex data (BPCD). 
% Modified based on Toll et al.(2020) and Hipp et al.(2013)
% BPCD = band passed complex data, chans (or verts) x timepoints, e.g., 64 chans x (180 sec * EEG.srate)
% RX = square correlation matrix for raw power envelopes, verts x verts or chans x chans, e.g., 5000 verts x 5000 verts
% RO = square correlation matric for orthog power envelopes, verts x verts or chans x chans, e.g., 5000 verts x 5000 verts
function [Rplain,Rortho] = OrthogonalPowCorr(BPCD)
   
nVerts = size(BPCD,1);
SS = single(BPCD);

SScda = conj(SS) ./ abs(SS); % SScda = Source Space conjugate divide by absolute;
PEX = SS .* conj(SS); % PEX = Power Envelopes Raw
PEX = log(PEX);
Rplain = single(corr(PEX',PEX'));

SScda = single(SScda);
PEX = single(PEX);
% RA is the orthogonal matrix;
RA = single(NaN(nVerts,nVerts));
waitbarHandle = waitbar(0,'Initializing','Name','Orthogonalized connectivity');
tic;
disp('Orthogonal power correlation calculation starts...');

for v=1:nVerts
  PEO = (imag(SS .* SScda(v,:) )).^2; % PEO = Power Envelopes Orthog
  PEO = log(PEO);
  PEXv = PEX(v,:);
  % correlate, there's a lot of overhead in matlab's corr function so the equation here is just the arithmetic for speed
  PEO = PEO - mean(PEO,2);
  PEXv = PEXv - mean(PEXv,2);
  RA(:,v) = sum(PEO.*PEXv,2) ./ ((sum(PEO.^2,2)).^.5 .* (sum(PEXv.^2,2)).^.5);
end

% A small value (tol) is added to the power envelopes to 
% prevent the natural logarithm resulting in a non-finite value (i.e., to prevent
% log(0) = -Inf). 
tol = single(1e-9); 

RA = gather(RA);
analysistime = toc;
disp(['The orthogonal power correlation analysis took ' num2str(analysistime) 's']);
delete(waitbarHandle);

% 8. Calculate the symmetric, corrected orthogonalized connectivity matrix % (Rortho). For region of interest intersection connectivity analyses, a
% symmetric connectivity matrix is required. Therefore, RorthoAB is averaged % with its transpose. Due to underestimation of correlation inherent to
% orthogonalization, a correction factor of 0.578499 is applied. See Hipp,
% J.F. et al., 2012. Large-scale cortical correlation structure of spontaneous % oscillatory activity. Nature Neuroscience, 15(6), pp.884?890. for the basis of % this.
Rortho = ((RA + RA') ./ 2) ./ 0.578499;
% 9. Limit the correlation values of Rplain and Rortho to |1 - tol|. This is % done to prevent non-finite values in the Fisher R to Z transform of the
% correlation coefficients (i.e., fisherz(1) = Inf).
Rplain(Rplain >= 1) = 1 - tol;
Rplain(Rplain <= -1) = -1 + tol; 
Rortho(Rortho >= 1) = 1 - tol; 
Rortho(Rortho <= -1) = -1 + tol;

% change the diagnal of Rplain to 0
Rplain(logical(eye(size(Rplain)))) = 0;

% r to z with two different ways that will get the same results;
Rortho = atanh(Rortho);
%Rorthoz = 0.5*log((1+Rortho)./(1-Rortho));
Rplain = atanh(Rplain);
%Rplainz = 0.5*log((1+Rplain)./(1-Rplain));

%% The correlation between the power envelopes of the orthogonalized signals can also be calculated with the following codes
% the codes below are a bit slower but more straight forward. They will get the same results;
%{ 
for v=1:nVerts
  % Seed is the vertex that SS is being orthogonalized with respect to.
% Dimensionality is 1 x t time points.
   seed = SS(v,:);

% Calculate PEseed by multiplication with its conjugate.
% Dimensionality is 1 x t time points.
  PEseed = seed .* conj(seed);  
% Calculate PEortho by applying Hipp's equation.
% Dimensionality is v vertices x t time points.
 PEortho = imag(SS .* ( conj(seed) ./ abs(seed) )).^2;
 % Calculate the correlation between the natural logarithm of PEseed 
 % and the natural logarithm of PEortho. This produces a column vector 
 % of v correlation coefficients per iteration.
 RA(:,v) = corr(log(PEseed' + tol),log(PEortho' + tol)); 
end
%}