function plotsession(data,multiunitdata,plottype)

switch plottype
%     case 'csd'
%         t = data.lfps.time;
%         ves_csd = data.csd.ves;
%         vis_csd = data.csd.vis;
%         nch = size(ves_csd,1);
%         % vestibular
%         z = 0.1*(2:nch+1);
%         [zz,tt] = meshgrid(z,t);
%         z_smooth = linspace(min(z),max(z),1000);
%         ves_csdsmooth = interp2(zz,tt,ves_csd',z_smooth,t,'spline');
%         figure; set(gcf,'Position',[100 100 400 700]);
%         imagesc(t,z_smooth,ves_csdsmooth');
%         hold on;
%         for i=1:nch
%             t2 = multiunitdata(i+1).ves.time;
%             r2 = mean(multiunitdata(i+1).ves.rate_pst);
%             r2 = r2 - mean(r2(t2>0 & t2<0.4));
%             r2 = r2/max(r2);
%             plot(t2,0.1*(i+1.5 - r2),'linewidth',2,'Color','k');
%         end
%         xlim([0 1.5]); colormap(brewermap([],'PuOr'));
%         % visual
%         figure; set(gcf,'Position',[500 100 400 700]);
%         z = 0.1*(2:nch+1);
%         [zz,tt] = meshgrid(z,t);
%         z_smooth = linspace(min(z),max(z),1000);
%         vis_csdsmooth = interp2(zz,tt,vis_csd',z_smooth,t);
%         imagesc(t,z_smooth,vis_csdsmooth');
%         hold on;
%         for i=1:nch
%             t2 = multiunitdata(i+1).vis.time;
%             r2 = mean(multiunitdata(i+1).vis.rate_pst);
%             r2 = r2 - mean(r2(t2>0 & t2<0.4));
%             r2 = r2/max(r2);
%             plot(t2,0.1*(i+1.5 - r2),'linewidth',2,'Color','k');
%         end
%         xlim([0 1.5]); colormap(brewermap([],'PuOr'));
%         % spatial correlation of lfp
%         figure; set(gcf,'Position',[900 400 400 300]);
%         hold on; 
%         dx = data.lfps.corr_distance.dx;
%         rves_mu = data.lfps.corr_distance.ves.r_mu;
%         rvis_mu = data.lfps.corr_distance.vis.r_mu;
%         plot(dx,rves_mu,'linewidth',3,'Color','r');
%         plot(dx,rvis_mu,'linewidth',3,'Color','g');
%         ylim([-1 1]);
    case 'csd'
        ves_csd = data.csd.ves;
        vis_csd = data.csd.vis;
        t = data.lfps.time;
        % vestibular
        nch = size(ves_csd,1);
        z = 0.1*(2:nch+1);
        [zz,tt] = meshgrid(z,t);
        z_smooth = linspace(min(z),max(z),1000);
        ves_csdsmooth = interp2(zz,tt,ves_csd',z_smooth,t,'spline');
        [ves_z,ves_t] = getsink(ves_csdsmooth');
        for i=1:nch
            t2 = multiunitdata(i+1).ves.time;
            r2(i,:) = mean(multiunitdata(i+1).ves.rate_pst);
            r2(i,:) = r2(i,:) - mean(r2(i,t2>0 & t2<0.4));
            r2(i,:) = r2(i,:)/max(r2(i,:));
        end
        [zz,tt2] = meshgrid(z,t2);
        ves_muasmooth = interp2(zz,tt2,r2',z_smooth,t,'spline');
        % visual
        nch = size(ves_csd,1);
        z = 0.1*(2:nch+1);
        [zz,tt] = meshgrid(z,t);
        z_smooth = linspace(min(z),max(z),1000);
        vis_csdsmooth = interp2(zz,tt,vis_csd',z_smooth,t,'spline');
        [vis_z,vis_t] = getsink(vis_csdsmooth');
        clear r2;
        for i=1:nch
            t2 = multiunitdata(i+1).vis.time;
            r2(i,:) = mean(multiunitdata(i+1).vis.rate_pst);
            r2(i,:) = r2(i,:) - mean(r2(i,t2>0 & t2<0.4));
            r2(i,:) = r2(i,:)/max(r2(i,:));
        end
        [zz,tt2] = meshgrid(z,t2);
        vis_muasmooth = interp2(zz,tt2,r2',z_smooth,t,'spline');
        load csd.mat;
        sessions(end+1).ves_csd = ves_csdsmooth;
        sessions(end).ves_mua = ves_muasmooth;
        sessions(end).ves_zsink= ves_z;
        sessions(end).ves_tsink= ves_t;
        sessions(end).vis_csd = vis_csdsmooth;
        sessions(end).vis_mua = vis_muasmooth;
        sessions(end).vis_zsink = vis_z;
        sessions(end).vis_tsink = vis_t;
        
        sessions(end).t_csd = t;
        sessions(end).t_mua = t2;
        sessions(end).z = z_smooth;
        save('csd.mat','sessions');
    case 'response_lfp'
    case 'response_multiunit'
    case 'response_singleunit'
    case 'corr_lfp'
end