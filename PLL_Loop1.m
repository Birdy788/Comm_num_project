function [ xI ] = PLL_Loop1(TED, L, mfIn, mfOut, K1, K2,const, Ksym)
%% Symbol Timing Loop
% ---------------------
%
% Implements closed-loop symbol timing recovery with a configurable timing
% error detector (TED) and linear interpolator. The feedback
% control loop uses a proportional-plus-integrator (PI) controller and a
% modulo-1 counter to control the interpolator. Meanwhile, the TED can be
% configured from 3 alternative implementations:
%     - Early-late TED (ELTED);
%     - Zero-crossing TED (ZCTED);
%     - Gardner TED (GTED);

% The interpolator used in this code is the Linear interpolator

% note that all TED schemes except the GTED compute the symbol
% timing error using symbol decisions. Hence, the reference constellation
% and scaling factor must be provided through input arguments 'const' and
% 'Ksym'. Such TED schemes could also leverage prior knowledge and compute
% the timing error using known symbols instead of decisions. However, this
% function does not offer the data-aided alternative. Instead, it only
% implements the decision-directed flavor of each TED. Meanwhile, the
% GTED is the only scheme purely based on the raw input samples, which is
% not decision-directed nor data-aided. Hence, the 'const' and 'Ksym'
% arguments are irrelevant when using the GTED.
%
% Input Arguments:
% TED     -> TED scheme ('MLTED', 'ELTED', 'ZCTED', 'GTED', or 'MMTED').
% L       -> Oversampling factor. ( It should always be equal to 1)
% mfIn    -> MF input sequence sampled at L samples/symbol.
% mfOut   -> MF output sequence sampled at L samples/symbol.
% K1      -> PI controller's proportional gain.
% K2      -> PI controller's integrator gain.
% const   -> Symbol constellation.
% Ksym    -> Symbol scaling factor to be undone prior to slicing. ( Only
% used of ZCTED & ELTED).

midpointOffset = ceil(L / 2);
muOffset = midpointOffset - L/2; % 0.5 if L is odd, 0 if L is even

% Modulation order
M = numel(const);

% Make sure the input vectors are column vectors
if (size(mfIn, 1) == 1)
    mfIn = mfIn(:);
end
if (size(mfOut, 1) == 1)
    mfOut = mfOut(:);
end

% Create an alias for the input vector to be used. As explained above, use
% the MF input with the polyphase interpolator and the MF output otherwise.

inVec = mfOut;

%% Timing Recovery Loop
polyMf = []; % dummy variable (used for the linear interpolator)
b_mtx = []; % dummy variable (used for the linear interpolator)

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

if (strcmp(TED, 'ELTED'))
    n_end = nSamples - L;
else
    n_end = nSamples;
end

for n = n_start:n_end
    if strobe == 1
        % Interpolation
        xI(k) = interpolate(inVec, m_k, mu(k), b_mtx, polyMf);
         % Timing Error Detector:
       % a_hat_k = Ksym * slice(xI(k) / Ksym, M); % Data Symbol Estimate
        switch (TED)
             case 'ELTED' % Early-late TED
                % Early and late interpolants
                early_idx = m_k + midpointOffset;
                late_idx = m_k - midpointOffset;
                early_mu = mu(k) - muOffset;
                late_mu = mu(k) + muOffset;
                x_early = interpolate( inVec, early_idx, ...
                    early_mu, b_mtx, polyMf);
                x_late = interpolate( inVec, late_idx, ...
                    late_mu, b_mtx, polyMf);
                % Decision-directed version of (8.99), i.e., (8.34)
                % adapted to complex symbols:
                e(n) = real(a_hat_k) * (real(x_early) - real(x_late)) + ...
                    imag(a_hat_k) * (imag(x_early) - imag(x_late));
            case 'ZCTED' % Zero-crossing TED
                % Estimate of the previous data symbol
                a_hat_prev = Ksym * slice(last_xI / Ksym, M);

                % Zero-crossing interpolant
                zc_idx = m_k - midpointOffset;
                zc_mu = mu(k) + muOffset;
                x_zc = interpolate( inVec, zc_idx, zc_mu, ...
                    b_mtx, polyMf);

                % Decision-directed version of (8.100), i.e., (8.37)
                % adapted to complex symbols:
                e(n) = real(x_zc) * ...
                    (real(a_hat_prev) - real(a_hat_k)) + ...
                    imag(x_zc) * (imag(a_hat_prev) - imag(a_hat_k));
            case 'GTED' % Gardner TED
                % Zero-crossing interpolant, same as used by the ZCTED
                zc_idx = m_k - midpointOffset;
                zc_mu = mu(k) + muOffset;
                x_zc = interpolate( inVec, zc_idx, zc_mu, ...
                    b_mtx, polyMf);

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
function [xI] = interpolate( x, m_k, mu, b_mtx, poly_f)
% [xI] = interpolate(method, x, m_k, mu, b_mtx, poly_h) returns the
% interpolant xI obtained from the vector of samples x.
%
% Args:
%     x      -> Vector of samples based on which the interpolant shall be
%               computed, including the basepoint and surrounding samples.
%     m_k    -> Basepoint index, the index preceding the interpolant.
%     mu     -> Estimated fractional interval between the basepoint index
%               and the desired interpolant instant.
%     b_mtx  -> Matrix with the coefficients for the polynomial
%               interpolator used with method=2 or method=3.
%     poly_f -> Polyphase filter bank that should process the input samples
%               when using the polyphase interpolator (method=0).

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
    xI = mu * x(m_k + 1) + (1 - mu) * x(m_k); % Linear Interpolator (See Eq. 8.61)
end

%This function is only used for the Early-late TED (ELTED)or the Zero-crossing TED (ZCTED). In our case we'll be using a GTED so this fucntion is not important for us but i just leave it in the code in case we wanna use 
%
%% Function to map Rx symbols into constellation points
% function [z] = slice(y, M)
% if (isreal(y))
%     % Move the real part of input signal; scale appropriately and round the
%     % values to get ideal constellation index
%     z_index = round( ((real(y) + (M-1)) ./ 2) );
%     % clip the values that are outside the valid range
%     z_index(z_index <= -1) = 0;
%     z_index(z_index > (M-1)) = M-1;
%     % Regenerate Symbol (slice)
%     z = z_index*2 - (M-1);
% else
%     M_bar = sqrt(M);
%     % Move the real part of input signal; scale appropriately and round the
%     % values to get ideal constellation index
%     z_index_re = round( ((real(y) + (M_bar - 1)) ./ 2) );
%     % Move the imaginary part of input signal; scale appropriately and
%     % round the values to get ideal constellation index
%     z_index_im = round( ((imag(y) + (M_bar - 1)) ./ 2) );
% 
%     % clip the values that are outside the valid range
%     z_index_re(z_index_re <= -1)       = 0;
%     z_index_re(z_index_re > (M_bar-1)) = M_bar-1;
%     z_index_im(z_index_im <= -1)       = 0;
%     z_index_im(z_index_im > (M_bar-1)) = M_bar-1;
% 
%     % Regenerate Symbol (slice)
%     z = (z_index_re*2 - (M_bar-1)) + 1j*(z_index_im*2 - (M_bar-1));
% end
% end
end

