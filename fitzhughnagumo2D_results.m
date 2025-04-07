clear all
close all

methods = ["ERK2ds","ERK2","Lawson2","IMEX2_LU","IMEX2_Sylvester","IMEX2_Sylvester_eig","ode23tb","ode23"];

R = cell(1,length(methods));
T = R;
for ii = 1:length(methods)
    addpath("files/results/fitzhughnagumo2Dresults/"+methods(ii))
    files = dir("files/results/fitzhughnagumo2Dresults/"+methods(ii));
    files = {files.name};
    files = files(3:end);
    for jj = 1:length(files)
        results = load(files{jj});
        R{ii}(jj,1) = results.par;
        R{ii}(jj,2) = results.time;
        R{ii}(jj,3) = results.err;
        if ~strcmp(methods(ii),"ode23tb") && ~strcmp(methods(ii),"ode23")
            if jj ~= 1
                R{ii}(jj,4) = -log(R{ii}(jj,3)/R{ii}(jj-1,3))/log(R{ii}(jj,1)/R{ii}(jj-1,1));
            else
                R{ii}(jj,4) = NaN;
            end
        end
    end
    if ~strcmp(methods(ii),"ode23tb") && ~strcmp(methods(ii),"ode23")
        T{ii} = table(char(compose('%d',R{ii}(:,1))),char(compose('%.2f',R{ii}(:,2))),...
            char(compose('%.3e',R{ii}(:,3))),char(compose('%.2f',R{ii}(:,4))),...
            'VariableNames',{'m','time (s)','error','order'});
    else
        [~,id] = sort(R{ii}(:,1),'descend');
        R{ii} = R{ii}(id,:);
        T{ii} = table(char(compose('%.0e',R{ii}(:,1))),char(compose('%.2f',R{ii}(:,2))),...
            char(compose('%.3e',R{ii}(:,3))),'VariableNames',{'tolerance','time (s)','error'});
    end
    T{ii} = table(T{ii},'VariableNames',{char(methods(ii))});
    disp(T{ii})
end
figure
markers = {'o-','+-','*-','>-','<-','-v','x-','^-'};
hold on 
for ii = 1:length(methods)
    loglog(R{ii}(:,2),R{ii}(:,3),markers{ii},'LineWidth',2)
end
hold off
set(gca,'XScale','log','YScale','log');
legend(cellstr(strrep(methods,'_','\_')),'Location','bestoutside')
xlim([1e1-5,1e4-5])
ylim([1e-3,1e-1])