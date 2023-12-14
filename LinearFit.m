function [fitresult, gof] = LinearFit(a, b, w)
%CREATEFIT8(A,B,W)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : a
%      Y Output: b
%      Weights : w
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 05-Apr-2019 18:56:03


%% Fit: 'untitled fit 1'.
[xData, yData, weights] = prepareCurveData( a, b, w );

% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );

% ft = fittype( 'smoothingspline' );
% opts = fitoptions( 'Method', 'SmoothingSpline' );

opts.Weights = weights;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'b vs. a with w', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel a
% ylabel b
% grid on


