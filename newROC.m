
function CP = newROC(Preferred, NonPref,CountorZscore) % 1 means spike count, 0 means z-score


Prefn = numel(Preferred);  % correct trials
NonPrefn = numel(NonPref); % incorrect trials

low = min([min(Preferred) min(NonPref)])-1; 
high = max([max(Preferred) max(NonPref)])+1;

if 1 == CountorZscore
    z = low:high;
    N = numel(z);
else
    N = 500;
    z=linspace(low*1.1,high*1.1,N); 
end

for i = 1:N
  fa(N-i+1) = sum(NonPref > z(i)); 
  hit(N-i+1) = sum(Preferred > z(i));
end


fa = fa/NonPrefn;

hit = hit/Prefn;

%figure; plot(fa,hit);
CP = trapz(fa,hit);



