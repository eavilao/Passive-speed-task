% plot example 1DAzi neurons
count = 0;
for i = [9 14 28 33 41 64 67 72 92 95 100 108 118 119 130 131 132 142 143 146 148 150 162 166 179 181 182 184 187 190 194 195 204 211 214 216 217 243 246 251 265 267 268 291 299 318 322 323 325 327 342 347 355 356];
    experiments(3).plotunit('psth','singleunits',i);
    count = count+1;
    hFig = gcf;
    print('1DAzi_Units','-append', '-dpsc2', '-bestfit'); % Print all
    % print('1DAzi_Units', '-dpdf', '-bestfit') % Print only one
    
    
    fprintf(['        Saved ' num2str(count)   '  ---> 1DAzi_Plots.pdf       \n'  ]);
    close all
end
