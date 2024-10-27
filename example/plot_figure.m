% This file conducts pointwise and joint inference.

clear
close all
clc

addpath(genpath(fullfile('..')))

%% PARAMETERS
% joint inference type
% 'TS': time-simultenous / 'MV-TS': multivariate, time-simultenous
inference_type = 'MV-TS';

% joint inference central tendency
% 'bayes' / 'mean' / 'median'
cen_ten = 'median';

% crediblity level for pointwise and joint inference
cred = [0.68, 0.9];

% horizons to conduct inference
horizon_inference = 20;

%% ANALYSIS

% load data
load('estimation.mat')

% scaler for the multuplier response
y_g = mean(exp(data.log_rgdp_cap/100)./exp(data.log_rgov_exp_cap/100));

% put draws in appropriate format
irf_draws(:,1:3,:) = irf(1:horizon_inference,1:3,:);

% rescale irfs
irf_draws = irf_draws * std(data{:,1}/100)  /  median(irf_draws(1,1,:));

% calculate multiplier
multiplier_draws = y_g * cumsum(irf_draws(:,2,:),1) ./ cumsum(irf_draws(:,1,:),1);

% add multiplier to irf draws
all_draws = cat(2, multiplier_draws, irf_draws);

% drop government revenue irfs
all_draws(:,end,:) = [];

% pointwise inference
pointwise = pointwise_inference(all_draws, 'Credibility', cred);

% joint inference
% Error bands: Sidak
joint_sidak = joint_inference(all_draws, ...
    'Method', 'Sidak', 'Credibility', cred, ...
    'InferenceType', inference_type);

% Error bands: Bonferroni
joint_bonferroni = joint_inference(all_draws, ...
    'Method', 'Bonferroni', 'Credibility', cred, ...
    'CentralTendency', cen_ten, ...
    'InferenceType', inference_type);

% Error bands: sup-t (PQO)
joint_supt = joint_inference(all_draws, ...
    'Method', 'sup-t', 'Credibility', cred, ...
    'CentralTendency', cen_ten, ...
    'Calibration', 'PQO', 'InferenceType', inference_type);

% Error bands: min-max'
joint_min_max = joint_inference(all_draws, ...
    'Method', 'min-max', 'Credibility', cred, ...
    'LossFunction', 'absolute', 'CentralTendency', cen_ten, ...
    'InferenceType', inference_type);

% Error bands: min-max (LQO)
joint_min_max_lqo = joint_inference(all_draws, ...
    'Method', 'min-max', 'Credibility', cred, ...
    'LossFunction', 'absolute', 'CentralTendency', cen_ten, ...
    'Calibration', 'LQO', 'InferenceType', inference_type);

% Error bands: min-max (BDR)
joint_min_max_bdr = joint_inference(all_draws, ...
    'Method', 'min-max', 'Credibility', cred, ...
    'LossFunction', 'absolute', 'CentralTendency', cen_ten, ...
    'Calibration', 'BDR', 'InferenceType', inference_type);

% auxilary varibles
var_names = {'multiplier', 'government expenditure', 'gross domestic product'};
y_labels = {'units', 'percentage', 'percentage'};
y_lims = [0, 1.6; -1, 2; -0.2, 0.4];
y_ticks{1} = y_lims(1,1):0.2:y_lims(1,2);
y_ticks{2} = y_lims(2,1):0.5:y_lims(2,2);
y_ticks{3} = y_lims(3,1):0.10:y_lims(3,2);

my_blue     = [0.00, 0.45, 0.75];
my_purple   = [0.36, 0.16, 0.60];
my_orange   = [0.94, 0.61, 0.40];
my_brown    = [0.49, 0.25, 0.00];
my_red      = [0.81, 0.02, 0.02];
my_lila     = [1.00, 0.00, 1.00];
my_beige    = [0.81, 0.68, 0.68];
my_green    = [0.78, 0.88, 0.17];

%% FIGURE: Joint Inference
fig1 = figure('Position', [1.00, 1.00, 650, 900]);
tl = tiledlayout(3, 2);
tl.Padding = 'loose';
% tiledlayout('horizontal', 'TileSpacing', 'compact')
for vv = 1:3
    nexttile(2*vv - 1)
    hold on
    plot_bonferroni = plot(squeeze(joint_bonferroni.error_band68(:,vv,:)),...
        'Color', my_lila, 'LineWidth', 2, 'LineStyle', '-.');
    plot_sidak = plot(squeeze(joint_sidak.error_band68(:,vv,:)),...
        'Color', my_green, 'LineWidth', 2, 'LineStyle', '--');
    plot_supt = plot(squeeze(joint_supt.error_band68(:,vv,:)),...
        'Color', my_red, 'LineWidth', 2, 'LineStyle', '--');
    plot_min_max = plot(squeeze(joint_min_max.error_band68(:,vv,:)),...
        'Color', my_brown, 'LineWidth', 2, 'LineStyle', '--');
    plot_min_max_lqo = plot(squeeze(joint_min_max_lqo.error_band68(:,vv,:)),...
        'Color', my_beige, 'LineWidth', 2, 'LineStyle', '--');
    plot_min_max_bdr = plot(squeeze(joint_min_max_bdr.error_band68(:,vv,:)),...
        'Color', my_orange, 'LineWidth', 2, 'LineStyle', '--');
    plot_pointwise = plot(squeeze(pointwise.error_band68(:,vv,:)),...
        'Color', my_purple, 'LineWidth', 2, 'LineStyle', '-');

    xlim([1, 20]);
    if vv == 1
        yline(1, 'LineWidth', 2)
    else
        yline(0, 'LineWidth', 2)
    end
    grid on

    xticks([1, 4:4:horizon_inference]);
    xlabel('quarters')
    ylim(y_lims(vv,:))
    yticks(y_ticks{vv});
    ylabel(y_labels{vv})
    if vv == 1
        title({'- 68% credible region -'; var_names{vv}}, 'FontWeight', 'normal')
    else
        title(var_names{vv}, 'FontWeight', 'normal')
    end

end

%% Fig 2
for vv = 1:3
    nexttile(2*vv)
    hold on
    plot_bonferroni = plot(squeeze(joint_bonferroni.error_band90(:,vv,:)),...
        'Color', my_lila, 'LineWidth', 2, 'LineStyle', '-.');
    plot_sidak = plot(squeeze(joint_sidak.error_band90(:,vv,:)),...
        'Color', my_green, 'LineWidth', 2, 'LineStyle', '--');
    plot_supt = plot(squeeze(joint_supt.error_band90(:,vv,:)),...
        'Color', my_red, 'LineWidth', 2, 'LineStyle', '--');
    plot_min_max = plot(squeeze(joint_min_max.error_band90(:,vv,:)),...
        'Color', my_brown, 'LineWidth', 2, 'LineStyle', '--');
    plot_min_max_lqo = plot(squeeze(joint_min_max_lqo.error_band90(:,vv,:)),...
        'Color', my_beige, 'LineWidth', 2, 'LineStyle', '--');
    plot_min_max_bdr = plot(squeeze(joint_min_max_bdr.error_band90(:,vv,:)),...
        'Color', my_orange, 'LineWidth', 2, 'LineStyle', '--');
    plot_pointwise = plot(squeeze(pointwise.error_band90(:,vv,:)),...
        'Color', my_purple, 'LineWidth', 2, 'LineStyle', '-');

    xlim([1, 20]);
    if vv == 1
        yline(1, 'LineWidth', 2)
    else
        yline(0, 'LineWidth', 2)
    end
    if vv == 3
        lgd = legend([plot_bonferroni(1), plot_sidak(1), plot_supt(1), ...
            plot_min_max(1), plot_min_max_lqo(1), plot_min_max_bdr(1),...
            plot_pointwise(1)], ...
            {'Bonferroni', 'Šidák', 'sup-t (PQO)', 'min-max', ...
            'min-max (LQO)', 'min-max (BDR)', 'pointwise'}, ...
            'NumColumns', 4, 'Box', 'off');
        lgd.Layout.Tile = 'south';
    end
    xticks([1, 4:4:horizon_inference]);
    xlabel('quarters')
    ylim(y_lims(vv,:))
    yticks(y_ticks{vv});
    ylabel(y_labels{vv})
    if vv == 1
        title({'- 90% credible region -'; var_names{vv}}, 'FontWeight', 'normal')
    else
        title(var_names{vv}, 'FontWeight', 'normal')
    end
    grid on
end
fontsize(fig1, scale=1.25)