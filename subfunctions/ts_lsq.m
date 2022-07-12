function [a,b,std_yx,hypoth,s_cf_ab,XYarray,r]=ts_lsq(x,y,alph,ISFIG,Ax,Leg,M,b_true)
%
% LEAST SQUARES STRAIGHT LINE REGRESSION
%
%            Y  =  ABSCISSA  +  SLOPE * X
%
%
% Syntax:  [ABSCISSA       ,... 
%           SLOPE          ,... 
%           STDEV          ,... 
%           HYPOTHESIS     ,... 
%           STD_CF_AB      ,... 
%           XYMATRIX       ,...
%           R              ]=mylsq(XDATA,YDATA,SIGNIFICANCE_LEVEL, ISFIG,Ax,Leg,XVECTOR,BTRUE)
%
% or just...
%
% Syntax:  [a , b ]=mylsq(x,y,alph,ISFIG,Ax,Leg)
% Syntax:  [a,b,s_yx,hypoth,s_cf_ab,XYarray,r]=mylsq(x,y,alph,ISFIG,Ax,Leg,M,btrue)
%
% INPUT:   XDATA     - The x-ordinates
%          YDATA     - The y-ordinates
%          SIGNIFICANCE_LEVEL - The significance level of the statistical test
%                               to be carried out. Example: 0.05 significance level.
%          ISFIG     - Use 1 for plotting on a figure. set to string 'silent'
%                      to omit any output to terminal. ('mute','quiet')
%          Ax        - If ISFIG>2 and Ax=axis handle then the plot will appear on that axis
%          Leg       - Legend (use 0)
%          XVECTOR   - A user defined monotonically increasing row-vector of X-values for the
%                      plotting and output of standard deviations and confidence limits.
%                      Ideally it could be XVECTOR = linspace(min(XDATA),max(XDATA),20)';
%          BTRUE     - The true slope of the population regression line. 
%                      Used for hypothesis testing
%
% OUTPUT:  ABSCISSA  - The abscissa of the least squares line
%          SLOPE     - The slope of the least squares line
%          STDEV     - The "standard error of estimate", that is, a measure
%                      of the "scatter in the vertical (y) direction of the
%                      observed points about the regression line".
%                      (Standard deviation of the regression line in the y-direction)
%          HYPOTHESIS- Vector of hypothesis testing values:
%                      1 = Degrees of freedom (n-2) (from input)
%                      2 = SIGNIFICANCE_LEVEL, the postulated significance level 
%                      3 = The corresponding critical value of the twosided t-distribution for the given confidence
%                          level and degrees of freedom.
%                          Also:   t_(alph/2,df-2) [in MATLAB: >> tinv( (1-alph)+alph/2 ,df-2 );]
%                      4 = BTRUE: Null Hypothesis: SLOPE == BTRUE, the postulated value of the population slope (from input)
%                      5 = The corresponding value of the t-statistic.
%                      6 = The sf.l. required for the rejection of the null hypothesis.
%                      7 = The corresponding N.H. rejection confidence percentage.
%          STD_CF_AB  - A 3x2 matrix:
%                       STD_CF_AB(1,1 or 2) = The sample standard error of the ABSCISSA or SLOPE, respectively.
%                       STD_CF_AB(2,1 or 2) = The upper 100*(1-alph)% confidence limit for the 
%                                             true abscissa or slope, respectively, of the true 
%                                             regression line.
%                       STD_CF_AB(3,1 or 2) = The corresponding lower confidence limits.
%
%          XYMATRIX  - A matrix with the following columns:
%                      Col. 1   - The user defined vector XVECTOR
%                      Col. 2   - The corresponding  Y = a+b*XVECTOR  values of the regression line
%                      Col. 3,1 - The  Y+1*STDEV  values
%                      Col. 3,2 - The  Y-1*STDEV  values
%                      Col. 4,1 - The  Y+2*STDEV  values
%                      Col. 4,2 - The  Y-2*STDEV  values
%                      Col. 5,1 - The 100*(1-alph)% (upper) confidence limits for the regression ordinate.
%                      Col. 5,2 - The 100*(1-alph)% (lower) confidence limits for the regression ordinate.
%                                 ("confidence limits for an ordinate to the true regression line")
%          R         - Sample correlation coefficient. R^2 gives the fraction of variance of YDATA
%                      that is accounted for by the regression of y on x (the regression line).
%                                 
%
% EXAMPLE: X=(1:50)';
%          Y=rand(50,1);
%          [a,b,s_yx,hypoth,std_cd_ab,XY,R]=mylsq(X,Y,0.05,1,1,0,(-7:57)');
%
%
% REFERENCE: CROWE ET AL. (1960): STATISTICS MANUAL
%
%  WRITTEN : Benedikt Halldorsson, Ph.D.-Candidate. 
%            Dept. of Civil, Structural and Environmental Engineering
%            State University of New York at Buffalo
%            E-mail: bh25@buffalo.edu
%
%  DATE    : 04/04/2004             VERSION:  1.0 
%----------------------------------------------------------------------
% MOD: 2016-01-12 tsonne: added silent option, commented out unused code
%

% The slope of the least squares straight line
b=(x'-mean(x))*(y-mean(y)) / ((x'-mean(x))*(x-mean(x)));
% The abscissa of the least squares straight line
a=mean(y)-b*mean(x);

n=length(x);
%
% 
s_yx2=(n-1)/(n-2)*(var(y)-b^2*var(x));

% 2.074 is t for alpha/2 and N-2 corresponding to P(t)=0.025
TALK = true;
if ischar(ISFIG)
    if any(strcmpi(ISFIG,{'no','silent','quiet','-q','mute'}))
        TALK = false; % no terminal output wanted
        ISFIG = false;
    end
else
    if ISFIG==1,
        fig;
    elseif ISFIG>1
        figure(gcf); axes(Ax); %subplot(ISFIG,1,ISFIG);
    end
end

if ISFIG
%h(1) = plot(x,y,'rd','markersize',6,'markerfacecolor','r');
h(1) = plot(x,y,'ko','markersize',6,'markerfacecolor','none','linewidth',1.5);
set(gca,'box','on');
hold on
end

if ~exist('M','var')
   %M=linspace(min(x),max(x),2*n)'; % tsonne: commented out
   M=linspace(min(x),max(x),25)'; % tsonne: subst.
end
if ~exist('b_true','var')
   b_true = 0;
end


XYarray(:,1,1) = M;
XYarray(:,2,1) = (a+b*M);
XYarray(:,3,1) = (a+b*M)+sqrt(s_yx2);
XYarray(:,3,2) = (a+b*M)-sqrt(s_yx2);
XYarray(:,4,1) = (a+b*M)+2*sqrt(s_yx2);
XYarray(:,4,2) = (a+b*M)-2*sqrt(s_yx2);


if ISFIG
h(2) = plot(M,XYarray(:,2,1),'b','linewidth',2);
h(3) = plot(M,XYarray(:,3,1),'k-.');
       plot(M,XYarray(:,3,2),'k-.');
h(4) = plot(M,XYarray(:,4,1),'k:');
       plot(M,XYarray(:,4,2),'k:');

ylabel('y','fontsize',12)
xlabel('x','fontsize',12)
end

std_yx=sqrt(s_yx2);
s_yx=sqrt(s_yx2);
s_x = sqrt( 1/(n-1) * (x-mean(x))'*(x-mean(x)) );
s_b = s_yx / ( s_x*sqrt(n-1) );
s_a = sqrt( s_yx^2 *(1/n + mean(x)^2/((x'-mean(x))*(x-mean(x)))) );

if TALK
disp(sprintf('--------------- R E S U L T S -------------------------------------------------------'));
disp(sprintf('*The sample regression line of y on x is '));
disp(sprintf('    Y = [ %1.4f +/- %1.6f] + X*[ %1.4f +/- %1.6f)]   (s_yx=%1.5f)',a,s_a,b,s_b,std_yx));
disp(sprintf('*s_yx is the standard error of the regression line for the Y-coordinate')); 
%disp(sprintf('*The dot-dashed and dotted lines show the regression line plus/minus one and two'));
%disp(sprintf('standard deviations (1*s_yx=%1.5f and 2*s_yx=%1.5f), respectively',std_yx,2*std_yx));
%disp(sprintf('*The dashed lines are the %2.0f confidence limits for the ordinate to the ',(1-alph)*100));
%disp(sprintf('true regression line'))
disp(sprintf('*Below is the result of testing the null hypothesis that the slope of the '));
disp(sprintf(' population regression line is %1.3f at the %2.0f%% significance level:',b_true,(1-alph)*100));
end

%==========================================================================
% TESTING IF THE SLOPE OF THE REGRESSION LINE IS ZERO OR NOT (Crowetal1960)
%==========================================================================
%
% THE SAMPLE STANDARD ERROR OF THE REGRESSION COEFFICIENT  a  IS (Bulmer1979,page 214)
%
               s_a = sqrt( s_yx^2 *(1/n + mean(x)^2/((x'-mean(x))*(x-mean(x)))) );
               s_cf_ab(1,1) = s_a;
%
% THE SAMPLE STANDARD ERROR OF THE REGRESSION COEFFICIENT  b  (THE SLOPE) IS:
%
               s_b = s_yx / ( s_x*sqrt(n-1) );                  %page 160. See also p.214 in Bulmer
               s_cf_ab(1,2) = s_b;
%
%--------------------------------------------------------------------------
%
% TO TEST THE NULL HYPOTHESIS THAT THE SLOPE b^ OF THE TRUE OR POPULATION REGRESSION LINE
% HAS THE STATED VALUE  b_true :
               if ~exist('b_true','var')
                  b_true = 0;
               end
% COMPUTE THE t-STATISTIC:
               t_test = (b - b_true) / s_b ;
%
% AND "REJECT THE NULL HYPOTHESIS AT THE SIGNIFICANCE LEVEL alph IF |t_test|
% EXCEEDS THE CRITICAL VALUE t_critical = t_(alph/2,n-2)" FROM THE t-DISTRIBUTION:
%
               %t_critical=my_ttable(n-2,alph,'t');     %my function
               t_critical=tinv((1-alph)+alph/2,n-2);   %built in MATLAB function
%
%--------------------------------------------------------------------------
%
% FINDING THE SIGNIFICANCE LEVEL AT WHICH THE NULL HYPOTHESIS CAN BE ACCEPTED
%
% IN THE SAME WAY I LOOKED UP THE t_critical CORRESPONDING TO alph 
% IN THE INVERSE STUDENT's t-CUMULATIVE DENSITY FUNCTION
% USING tinv(alph/2,n-2) I CAN,
% FOR THE GIVEN t_test VALUE, LOOK UP THE CORRESPONDING alph_test VALUE 
% OF THE STUDENT's t CUMULATIVE DENSITY FUNCTION (DISTRIBUTION) USING:
%
               alph_test_upper = tcdf(abs(t_test),n-2);          % F(t_test)=P(t|t<t_test): The area under t-cdf below t_test 
               alph_test = (1 - alph_test_upper)*2;         % The two sided area at the tails
               t_test_conflevel = 100*(1-alph_test);
%
% WHERE t_test_conflevel IS THE "REJECTION CONFIDENCE" OF THE NULL HYPOTHESIS
%
% For output:
               hypoth(1,1) = n-2;
               hypoth(2,1) = alph;
               hypoth(3,1) = t_critical;
               hypoth(4,1) = b_true;
               hypoth(5,1) = t_test;
               hypoth(6,1) = alph_test;
               hypoth(7,1) = t_test_conflevel;
%
% Printing out:               
%               
if TALK
if abs(t_test)>t_critical
      disp(sprintf('  |t_test| > t_(%1.2f/2,%d-2) <=> |%1.3f| > %1.3f is valid.  ',alph,n,abs(t_test),t_critical));
      disp(sprintf('  => Reject the null hypothesis that slope is %1.2f at the %1.0f%% significance level.',b_true,alph*100));
      disp(sprintf('  More specifically: Reject the null hypothesis that slope is %1.2f at the %1.0f%% significance level.',b_true,alph_test*100));
%      disp(sprintf('  i.e. There is %1.0f%% risk of error that I accept the null hypothesis when it is not true.',alph_test*100));
      disp(sprintf('  In other words: Rejection confidence of slope being %1.2f is %1.0f%%.',b_true,(1-alph_test)*100));
   else
      disp(sprintf('  |t_test| > t_(%1.2f/2,%d-2) <=> |%1.3f| > %1.3f is NOT valid.  ',alph,n,abs(t_test),t_critical));
      disp(sprintf('  => CANNOT Reject the null hypothesis that slope is %1.2f at the %1.0f%% significance level.',b_true,100*alph));
      disp(sprintf('  More specifically: Reject the null hypothesis that slope is %1.2f at the %1.0f%% significance level.',b_true,alph_test*100));
%      disp(sprintf('  i.e. There is %1.0f%% risk of error that I reject the null hypothesis when it is actually true.',alph_test*100));
      disp(sprintf('  In other words: Rejection confidence of slope being %1.2f is %1.0f%%.',b_true,(1-alph_test)*100));
end
end
%
%
%--------------------------------------------------------------------------
%
% A 100*(1-alph)% CONFIDENCE INTERVAL FOR THE TRUE INTERCEPT a_true IS:
%
               s_cf_ab(2,1) = a + t_critical*s_a;
               s_cf_ab(3,1) = a - t_critical*s_a;
%
% A 100*(1-alph)% CONFIDENCE INTERVAL FOR THE TRUE SLOPE b_true IS:
%
               s_cf_ab(2,2) = b + t_critical*s_b;
               s_cf_ab(3,2) = b - t_critical*s_b;
%
%--------------------------------------------------------------------------
%
% FINDING THE CONFIDENCE LIMITS FOR THE ORDINATE TO THE TRUE REGRESSION LINE (CROWEETAL1960,6.1.4e)
%
% WE CAN SHOW HOW GOOD OUR ESTIMATE Y' (AT ANY GIVEN X') OF THE TRUE MEAN ORDINATE IS
% BY CALCULATING THE 100*(1-alph)% CONFIDENCE LIMITS
%

               alph_cf_lim_y_true_line = s_yx * sqrt( 1/n+(M-mean(x)).^2/(n-1)/var(x) );  %page 162
%
               XYarray(:,5,1) = (a+b*M) + t_critical * alph_cf_lim_y_true_line ;
               XYarray(:,5,2) = (a+b*M) - t_critical * alph_cf_lim_y_true_line ;
%
% Plotting
if ISFIG
h(5) = plot(M,XYarray(:,5,1) ,'--','markersize',2);
       plot(M,XYarray(:,5,2) ,'--','markersize',2);
end
%
%--------------------------------------------------------------------------
%
% CALCULATING THE SAMPLE CORRELATION COEFFICIENT
%
% THE BEST ESTIMATE OF THE POPULATION CORRELATION COEFFICIENT rho (WHICH SQUARED GIVES
% THE FRACTION OF THE POPULATION y VARIANCE ACCOUNTED FOR BY THE REGRESSION ON x)
% IS THE SAMPLE CORRELATION COEFFICIENT r
%
                %r = (x'-mean(x))*(y-mean(y)) / (n-1)/std(x)/std(y);
                r = b * std(x)/std(y);                      %page 158
%
% THE SQUARE ROOT OF THIS FRACTION GIVES THE FRACTION OF VARIANCE OF y THAT IS 
% ACCOUNTED FOR BY THE REGRESSION LINE.
%

if TALK
disp(sprintf('-------------------------------------------------------------------------------------'));
end


return




if Leg
AX=legend(h,'\eta',...
         sprintf('y = %5.3f %+6.4f*M_w',a,b),...
         sprintf('\\pm \\sigma_y = %5.3f',sqrt(s_yx2)),...
         sprintf('\\pm 2\\sigma_y = %5.3f',2*sqrt(s_yx2)),...
         3);
%         '95% conf. interval for y-ordinate',...
%         '95% conf. interval for regression line',...
LEG = findobj(AX,'type','text'); set(LEG,'FontSize',8)
end

