close all
clearvars

% read data
% FtStripe = [10 20 40 60 80];  
% FtGrating = [.5 1 5 10 20 40 60 80 100 120];   
FtGratingx2 = [5.625 11.25 56.25 112.5 225 450 675 900 1125 1350]; FtGratingx2 = FtGratingx2 * 0.0889;
FtGratingx4 = [8.4375 16.875 84.375 168.75 337.5 675 1012.5 1350 1687.5 2025]; FtGratingx4 = FtGratingx4 * 0.0444;
FtGratingx8 = [14.0625 28.125 140.625 281.25 562.5 1125 1687.5 2250 2812.5 3375]; FtGratingx8 = FtGratingx8 * 0.0222;
% for single stripe
VsmallStripe = [112.5 225 450 675 900];
VmediumStripe = [168.75 337.5 675 1012.5 1350];
VlargeStripe = [281.25 562.5 1125 1687.5 2250];
% for gratings.25 11
VsmallGrating = [5.625 11.25 56.25 112.5 225 450 675 900 1125 1350];
VmediumGrating = [8.4375 16.875 84.375 168.75 337.5 675 1012.5 1350 1687.5 2025];
VlargeGrating = [14.0625 28.125 140.625 281.25 562.5 1125 1687.5 2250 2812.5 3375];

% plot
figure; hold on; box on
xlabel('temporal frequency (Hz)')
ylabel('Velocity (°/s)')
set(gca,'XTick',[0 20 40 60 80 100 120],'XTickLabels',{'0','20','40','60','80','100','120'})
axis square
xlim([0 130])
ylim([0 3500])
% plot([10 10],[0 3500],'--','Color',[.4 .4 .4 .4],'LineWidth',1.2)
% plot([80 80],[0 3500],'--','Color',[.4 .4 .4 .4],'LineWidth',1.2)
% patch([10 10 80 80],[0 3500 3500 0],[0.9 0.9 0.9],'EdgeColor','none','FaceAlpha',.5)
% gratings
plot(FtGratingx2,VsmallGrating,'-o','Color',[.2 .2 .2],'LineWidth',1)
plot(FtGratingx4,VmediumGrating,'-o','Color',[.2 .2 .2],'LineWidth',1.5)
plot(FtGratingx8,VlargeGrating,'-o','Color',[.2 .2 .2],'LineWidth',2)
% single stripe
plot(FtGratingx2(4:8),VsmallStripe,'-o','Color',[.6 .6 .6],'LineWidth',1)
plot(FtGratingx4(4:8),VmediumStripe,'-o','Color',[.6 .6 .6],'LineWidth',1.5)
plot(FtGratingx8(4:8),VlargeStripe,'-o','Color',[.6 .6 .6],'LineWidth',2)

% save figures
savefig('Methods_Velocity_Frequency.fig')
print('Methods_Velocity_Frequency','-depsc','-r300','-tiff','-painters')
