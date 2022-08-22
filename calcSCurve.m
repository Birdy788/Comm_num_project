function [ normTauE, g ] = calcSCurve(TED, rollOff)
%% Compute the S-Curve of a given TED assuming data-aided operation
%
% [ normTauE, g ] = calcSCurve(TED, rollOff, rcDelay) returns the S-curve
% g(normTauE) of the chosen TED computed analytically based on closed-form
% expressions. It also returns the vector normTauE with the normalized time
% offset errors (within -0.5 to 0.5) where the S-curve is evaluated.
%    TED -> TED choice.
%    rollOff -> Rolloff factor.
%
% Note: this function assumes constant channel gain (e.g., after automatic
% gain control) and unitary average symbol energy.


%% Parameters and assumptions
% Oversampling factor
%
% Note the oversampling factor that is effectively used on a symbol timing
% recovery loop has nothing to do with the factor adopted here. Here, the
% only goal is to observe the S curve, and the S curve (a continuous time
% metric) should be independent of the oversampling factor. Hence, we use a
% large value to observe the S-curve with enough resolution.
L = 1e3;

% Assume constant gain and unitary average symbol energy
K  = 1;
Ex = 1;

%% TED Gain

% Zero-centered indices corresponding to one symbol interval
tau_e = -L/2 : L/2; % timing error in units of sample periods
normTauE = tau_e / L; % normalized timing error

%% S-Curve
if (strcmp(TED, 'GTED'))
     % Gardner - See Eq. (8.45)
        C = sin(pi * rollOff / 2) / (4 * pi * (1 - (rollOff^2 / 4)));
        g = (4 * K^2 * Ex) * C * sin(2 * pi * tau_e / L);
        % NOTE:
        % - Eq. (8.45) normalizes tau_e by Ts inside the sine argument.
        %   Here, in contrast, we normalize it by L (oversampling factor).
        % - In Eq. (8.45), the scaling constant (4 * K^2 * Ex) is also
        %   divided by Ts. Here, we don't divide it, not even by L. This
        %   has been confirmed empirically using the "simSCurve" script.
end

end