function [h_main, h_inset]=inset(main_handle, inset_handle,fs)

% The function plotting figure inside figure (main and inset) from 2 existing figures.
% inset_size is the fraction of inset-figure size, default value is 0.35
% The outputs are the axes-handles of both.
% 
% An examle can found in the file: inset_example.m
% 
% Moshe Lindner, August 2010 (C).

if nargin==2
    inset_size=0.35;
end

% inset_size=inset_size*.35;
figure
new_fig=gcf;
main_fig = findobj(main_handle,'Type','axes');
h_main = copyobj(main_fig,new_fig);

FS='FontSize';
set(h_main,'Position',get(main_fig,'Position'));inset_fig = findobj(inset_handle,'Type','axes');h_inset = copyobj(inset_fig,new_fig);
title('Noise intensity');


ylim([0 1.5]);xticks([0 10 20]);xticklabels([0 10 20]);%a=get(gca,'XTickLabel');set(gca,'XTickLabel',a,FS,fs);hold on;
yticks([0 1.27]);yticklabels([0 1.27]);%a=get(gca,'YTickLabel');set(gca,'YTickLabel',a,FS,fs);

ax=get(main_fig,'Position');
% set(h_inset,'Position', [.9*ax(1)+ax(3)-inset_size .5*ax(2)+ax(4)-inset_size inset_size inset_size])
set(h_inset,'Position', [0.5 0.5 0.35 0.35])

