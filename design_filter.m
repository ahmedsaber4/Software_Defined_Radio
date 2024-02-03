function [IR]=design_filter(Fs, Fc, N, windowType)
    % Fs: Sample frequency
    % Fc: Cutoff frequency
    % N: Filter order (odd for type I FIR)
    % windowType: String specifying the window function

    % Generate the ideal impulse response of the low-pass filter
    t = (-N/2:N/2) / Fs;
    ideal_response = 2 * Fc / Fs * sinc(2 * Fc * t);

    % Generate and plot the response for the specified window
    if strcmpi(windowType, 'hamming')
        window = hamming(N);
    elseif strcmpi(windowType, 'hann')
        window = hann(N);
    elseif strcmpi(windowType, 'blackman')
        window = blackman(N);
    elseif strcmpi(windowType, 'chebyshev')
        window=chebwin(N);
    else
        error('Invalid window type');
    end

    % Apply window to the ideal response
    filtered_response = ideal_response .* window;

    % Normalize the filter coefficients
    IR = filtered_response / sum(filtered_response);
end
