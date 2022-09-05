function [ xI ] = PLL_Loop(TED, L, mfOut, K1, K2)
%% Symbol Timing Loop
    %This Code is inspired as all the equations from the book : Digital
    %Communications a Discrete Time AProach : Michael Rice. 

% ---------------------
%
% Implements closed-loop symbol timing recovery with a Gardner TED (GTED) and linear interpolator.
% The feedback control loop uses a proportional-plus-integrator (PI) controller and a
% modulo-1 counter to control the interpolator. 

% The interpolator used in this code is the Linear interpolator

% Input Arguments:
% TED     -> TED scheme : GTED.
% L       -> Oversampling factor. ( 1)
% mfOut   -> MF(Matched Filter) output sequence sampled at L samples/symbol.
% K1      -> PI controller's proportional gain.
% K2      -> PI controller's integrator gain.

midpointOffset = ceil(L / 2);
muOffset = midpointOffset - L/2; % 0.5 if L is odd, 0 if L is even


% Make sure the input vectors are column vectors
if (size(mfOut, 1) == 1)
    mfOut = mfOut(:);
end

% Create an alias for the input vector to be used. We use the MF output in
% our case.

inVec = mfOut;

%% Timing Recovery Loop

% Constants
nSamples = length(inVec);
nSymbols = ceil(nSamples / L);

% Preallocate
xI = zeros(nSymbols, 1); % Output interpolants
mu = zeros(nSymbols, 1); % Fractional symbol timing offset estimate
v  = zeros(nSamples, 1); % PI output
e  = zeros(nSamples, 1); % Error detected by the TED

% Initialize
k      = 0; % interpolant/symbol index
strobe = 0; % strobe signal
cnt    = 1; % modulo-1 counter
vi     = 0; % PI filter integrator
n_start = 1;
last_xI = inVec(1);
n_end = nSamples;


for n = n_start:n_end
    if strobe == 1
        % Interpolation
        xI(k) = interpolate(inVec, m_k, mu(k));
         % Timing Error Detector:
        if(strcmp(TED,'GTED')) % Gardner TED
                % Zero-crossing interpolant, same as used by the ZCTED
                zc_idx = m_k - midpointOffset;
                zc_mu = mu(k) + muOffset;
                x_zc = interpolate( inVec, zc_idx, zc_mu);

                % Equation (8.101):
                e(n) = real(x_zc) * (real(last_xI) - real(xI(k))) ...
                    + imag(x_zc) * (imag(last_xI) - imag(xI(k)));
        end
        
        % Update the "last output interpolant" for the next strobe
        last_xI = xI(k);
    else
        % Make the error null on the iterations without a strobe. This is
        % equivalent to upsampling the TED output.
        e(n) = 0;
    end
    % Loop Filter
    vp   = K1 * e(n);        % Proportional
    vi   = vi + (K2 * e(n)); % Integral
    v(n) = vp + vi;          % PI Output
    % NOTE: since e(n)=0 when strobe=0, the PI output can be simplified to
    % "v(n) = vi" on iterations without a strobe. It is only when strobe=1
    % that "vp != 0" and that vi (integrator output) changes. Importantly,
    % note the counter step W below changes briefly to "1/L + vp + vi" when
    % strobe=1 and, then, changes back to "1/L + vi" when strobe=0. In the
    % meantime, "vi" remains constant until the next strobe.

    % Adjust the step used by the modulo-1 counter (see below Eq. 8.86)
    W = 1/L + v(n);

    % Check whether the counter will underflow on the next cycle, i.e.,
    % whenever "cnt < W". When that happens, the strobe signal must
    % indicate the underflow occurrence and trigger updates on:
    %
    % - The basepoint index: set to the index right **before** the
    %   underflow. When strobe=1, it means an underflow will occur on the
    %   **next** cycle. Hence, the index before the underflow is exactly
    %   the current index.
    % - The estimate of the fractional symbol timing offset: the estimate
    %   is based on the counter value **before** the underflow (i.e., on
    %   the current cycle) and the current counter step, according to
    %   equation (8.89).
    
    strobe = cnt < W;
    if (strobe)
        k = k + 1; % Update the interpolant Index
        m_k = n; % Basepoint index (the index **before** the underflow)
        mu(k) = cnt / W; % Equation (8.89)
    end

    % Next modulo-1 counter value:
    cnt = mod(cnt - W, 1);
end

% Trim the output vector
if (strobe) % ended on a strobe (before filling the k-th interpolant)
    xI = xI(1:k-1);
else
    xI = xI(1:k);
end

%% Interpolation
function [xI] = interpolate( x, m_k, mu)
% [xI] = interpolate(method, x, m_k, mu) returns the
% interpolant xI obtained from the vector of samples x.
%
% Args:
%     x      -> Vector of samples based on which the interpolant shall be
%               computed, including the basepoint and surrounding samples.
%     m_k    -> Basepoint index, the index preceding the interpolant.
%     mu     -> Estimated fractional interval between the basepoint index
%               and the desired interpolant instant.

% Adjust the basepoint if mu falls out of the nominal [0,1) range. This
% step is necessary only to support odd oversampling ratios, when a
% +-0.5 offset is added to the original mu estimate. In contrast, with
% an even oversampling ratio, mu is within [0,1) by definition.
    if (mu < 0)
        m_k = m_k - 1;
        mu = mu + 1;
    elseif (mu >= 1)
        m_k = m_k + 1;
        mu = mu - 1;
    end
    assert(mu >= 0 && mu < 1); 
    c0=mu;
    c1=1-mu;
    xI = c0 * x(m_k + 1) + c1 * x(m_k); % Linear Interpolator (See Eq. 8.61)
end

end

