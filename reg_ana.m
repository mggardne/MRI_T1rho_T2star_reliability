function [tc,amp,rss,npx,id,tcpx,amppx,rsspx,nps] = reg_ana(v,mask, ...
                             rsls,nrsls,irsl,time,ntime,fun,init,t0,opt)
%REG_ANA   Given a four-dimensional matrix of T1/T2 intensities from
%          a MRI image volume, a mask, and a vector of spin lock/echo
%          times, calculates the T1rho/T2* map for the volume.
%
%          TC = REG_ANA(V,MASK,RSLS,NRSLS,IRSL,TIME,NTIME,FUN,INIT)
%          Given a four-dimensional matrix of T1/T2 intensities from a
%          MRI image volume, V, where the first two dimensions are an
%          image, the third dimension are the slices and the fourth
%          dimension are the spin lock/echo times, three dimensional
%          logical masks with the first dimension being the image, the
%          second dimension being the superficial layer in the first
%          column and deep layer in the second column and the third
%          dimension being slices in a cell array of masks with the
%          first index to regions 1 and 2 and the second index to the
%          femur and tibia, MASK, a cell array with the slices within
%          each region, RSLS, the number of slices in each region,
%          NRSLS, an identifier for the compartment of each region,
%          IRSL, vector of spin lock/ echo times, TIME, the number of
%          spin lock/echo times, NTIME, a function handle, FUN, to
%          evaluate the monoexponential and its derivatives, an
%          initialization flag, INIT, to determine the initial starting
%          parameters, calculates T1rho/T2* for the two compartments
%          (lateral and medial), two bones (femur and tibia) and two
%          layers (superficial and deep) in array TC.
%
%          TC = REG_ANA(V,MASK,RSLS,NRSLS,IRSL,TIME,NTIME,FUN,INIT,T0)
%          If INIT is greater than zero, uses T0 as the initial default
%          time constant for the nonlinear exponential fit.  The
%          default value is 50.  The mean value of the maximum
%          intensities are used as the initial amplitude.  If INIT is
%          equal to zero, the exponential equation is linearized using
%          logarithms and least squares are used to solve for the
%          initial parameters for the nonlinear least squares.  If INIT
%          is less than zero, weighted least squares using the
%          intensities is used to solve for the starting parameters.
%
%          TC = REG_ANA(V,MASK,RSLS,NRSLS,IRSL,TIME,NTIME,FUN,INIT,T0,
%          OPT) a structure, OPT, with options for the curvefit solver.
%          See Matlab command "optimset" for a list of parameters.
%
%          [TC,AMP,RSS,NPX,ID] = REG_ANA(...) returns the amplitudes of
%          the exponential fits, AMP, the sum of squared residuals, RSS,
%          the number of pixels in the curvefits, NPX, and a three
%          column array with the first column being the compartment (0 -
%          lateral and 1 - medial), the second column being bone (0 -
%          femur and 1 - tibia), and the third column being layer (0 -
%          deep and 1 - superficial) that identifies the analyzes, ID. 
%
%          [TC,AMP,RSS,NPX,ID,TCPX,AMPPX,RSSPX,NPS] = REG_ANA(...)
%          returns the T1rho/T2* for each pixel in the two compartments
%          (lateral and medial), two bones (femur and tibia) and two
%          layers (superficial and deep) in a cell array, TCPX, the
%          exponential amplitudes for each pixel in a cell array, AMPPX,
%          the sum of squared residuals for each pixel, RSSPX, and the
%          number of pixels in each slice within the regions in a cell
%          array, NPS.
%
%          NOTES:  1.  The T1rho/T2* map is calculated as a mono-
%                  exponential as a function of spin lock/echo times.
%
%                  2.  The images volume data and spin lock/echo times
%                  are assumed to be in ascending order with increasing
%                  index into the fourth dimension of the image data V.
%
%                  3.  The function handle, FUN, is to the M-file
%                  exp_fun1.m which must be in the current directory or
%                  path.
%
%                  4.  T0 must be between 1 and 100.
%
%                  5.  The Matlab Optimization toolbox is required.
%                  The nonlinear least squares is performed by the
%                  Matlab function lsqcurvefit in the optimization
%                  toolbox.
%
%                  6.  The nonlinear least squares uses the Levenberg-
%                  Marquardt algorithm in the Matlab function
%                  lsqcurvefit.
%
%                  7.  The Matlab Parallel Computing toolbox is 
%                  required.  The Matlab parallel construct parfor is
%                  used to calculate the T1rho/T2* values in parallel.
%                  Use Matlab command parpool to control the number of 
%                  workers.
%
%                  8.  Region is NOT tibial compartment.  Variable IRSL
%                  indicates the tibial compartment of each region.  The
%                  relationship between region and compartment depends
%                  on the leg and the direction the slices were
%                  digitized.
%
%          26-Jan-2022 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if nargin<9
  error([' *** ERROR in reg_ana:  Nine input variables are', ...
         ' required!']);
end
%
if nargin<10
  t0 = 50;              % Initial estimate of time constant
end
%
% Check Initial Time Constant
%
if isempty(t0)||t0<1||t0>100
  t0 = 50;
end
%
% Default Optimization Parameters
%
if nargin<10
  opt = optimset('Display','off','TolFun',1e-8,'TolX',1e-8, ...
                 'MaxIter',2e+3,'Algorithm','levenberg-marquardt', ...
                 'Jacobian','on','UseParallel',true);
end
%
% Initialize Variables
%
% Columns in variable "id":
%   1 - Compartments:  0 - Lateral and 1 - Medial
%   2 - Bones:  0 - Femur and 1 - Tibia
%   3 - Layers:  0 - Deep and 1 - Superficial
%
tc = zeros(8,1);        % Time constant - T1rho/T2*
amp = zeros(8,1);       % Amplitude of exponential
rss = zeros(8,1);       % Sum of squared residuals
npx = zeros(8,1);       % Number of pixels in regional fit
id = zeros(8,3,1);      % Columns: 1 - Compartment, 2 - Bone and 3 - Layer
%
nps = cell(8,1);        % Number of pixels in each slice
tcpx = cell(8,1);       % Time constant - T1rho/T2* for each pixel
amppx = cell(8,1);      % Amplitude of exponential for each pixel
rsspx = cell(8,1);      % Sum of squared residuals for each pixel
%
% Loop through Regions
%
for kr = 1:2
%
   ikr = irsl(kr)-1;    % Compartment code (0 - lateral, 1 - medial)
   mskr = mask{kr};     % Mask for this region
   rsl = rsls{kr};      % Slices for this region
   nrsl = nrsls(kr);    % Number of slices in this region
%
% Loop through Bones
%
   for kb = 1:2
%
      ikb = kb-1;       % Bone code identifier (0 - femur, 1 - tibia)
      msk = mskr{kb};   % Mask for this bone
      rimgs = cell(nrsl,2);            % # of slices, # of layers
%
% Loop through Slices
%
      np = zeros(2,nrsl);
%
      for ks = 1:nrsl
%
         slk = rsl(ks);                % Slice number
%
         rimgt = cell(ntime,2);        % # of spin lock/echo times, # of layers
%
% Loop through Spin Lock/Echo Times
%
         for kt = 1:ntime
%
            rimg = squeeze(v(:,:,slk,kt));  % T1/T2 data for slice and spin lock/echo time
%
% Loop through Layers
%
            for kl = 1:2
               rimgt{kt,kl} = rimg(msk(:,kl,ks))';
               np(kl,ks) = sum(msk(:,kl,ks));
            end
%
         end            % End of kt loop - Spin lock/echo times
%
         for kl = 1:2
            rimgs{ks,kl} = cat(1,rimgt{:,kl});   % Combine spin lock/echo times
         end
         
      end               % End of ks loop - slices loop
%
      clear rimg rimgt;
%
% Calculate T1rho/T2* for Each Layer
%
      for kl = 1:2
%
         idx = 4*kr+2*kb+kl-6;         % Index to variables
         nps{idx} = np(kl,:);          % Number of pixels in slice
%
         rimgl = cat(2,rimgs{:,kl});
         npx(idx) = size(rimgl,2);
         r = rimgl(:);
%
         amp0 = mean(rimgl(1,:));
         rp0c = [amp0; t0];            % Constant initial parameters
%
         tl = repmat(time,npx(idx),1);
%
         if init==0     % Use linear least squares for initial constants
           xdat = [tl ones(ntime*npx(idx),1)];
           lr = log(r);
           rpl = xdat\lr;
           t0l = -1/rpl(1);
           if t0l>1&&t0l<100
             rp0 = [exp(rpl(2)); t0l];
             else
                rp0 = rp0c;
           end
         elseif init<0  % Use weighted least squares for initial constants
           xdat = [tl.*r r];
           lr = log(r);
           rpl = xdat\(lr.*r);
           t0l = -1/rpl(1);
           if t0l>1&&t0l<100
             rp0 = [exp(rpl(2)); t0l];
             else
                rp0 = rp0c;
           end
           else                % Use constant initial parameters
             rp0 = rp0c;
         end
%
% Nonlinear Least Squares Exponential Fit to Get T1rho/T2* Value
%
         [rp,err] = lsqcurvefit(fun,rp0,tl,r,[],[],opt);
%
         tc(idx) = rp(2);              % Time constant - T1rho/T2* 
         amp(idx) = rp(1);             % Amplitude of exponential
         rss(idx) = err;               % Residual sum of squared
         id(idx,:) = [ikr ikb 2-kl];   % Compartment, bone and layer IDs
%
% Calculate T1rho/T2* for Each Pixel in Each Layer
%
         if nargout>5
           tcp = zeros(npx(idx),1);
           ampp = tcp;
           rssp = tcp;
%
% Get Initial Parameters
%
           if init==0        % Use linear least squares for initial constants
             xdat = [time ones(ntime,1)];
             rpa = xdat\log(rimgl);
             t0a = -1./rpa(1,:);
             ampa = exp(rpa(2,:));
             rpa = [ampa; t0a];
             idv = t0a<1|t0a>100;
             if any(idv)
               rpa(:,idv) = repmat(rp0c,1,sum(idv));
             end
           elseif init<0     % Use weighted least squares for initial constants
             rpa = repmat(rp0c,1,npx(idx));
             for p = 1:npx(idx)
                rimgp = rimgl(:,p);
                xdat = [time.*rimgp rimgp];
                lr = log(rimgp);
                rpl = xdat\(lr.*rimgp);
                t0l = -1/rpl(1);
                if t0l>=1&&t0l<=100
                  rpa(:,p) = [exp(rpl(2)); t0l];
                end
             end
           else              % Use constant initial parameters
             rpa = repmat(rp0c,1,npx(idx));
           end
%
% Nonlinear Least Squares Exponential Fit to Get T1rho/T2* Pixel Values
%
           parfor p = 1:npx(idx)
%      
              rimgp = rimgl(:,p);
%      
              [rp,err] = lsqcurvefit(fun,rpa(:,p),time,rimgp,[],[],opt);
              tcp(p) = rp(2);
              ampp(p) = rp(1);
              rssp(p) = err;
%
           end
%
% Save Results to Outputs
%
           tcpx{idx} = tcp;
           amppx{idx} = ampp;
           rsspx{idx} = rssp;
%
         end
%
      end               % End of kl loop - layers loop
   end                  % End of kb loop - bones loop
end                     % End of kr loop - regions loop
%
return