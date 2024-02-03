clc; 
clear; 
close all;
% Number of NB-IoT devices
optimalfilter=cell2mat(struct2cell(load('3_Equiripple.mat')));
numDevices = 6;
IF_frequency_spacing = 200e3 ;  % Frequency spacing between NB-IoT carriers in Hz
IF_center_frequency = 0;  % Start at 0 Hz  % Center frequency for the composite IF signal in Hz
IF_carrier_freqs = IF_center_frequency + (0:numDevices-1) * IF_frequency_spacing;
% Display the calculated IF frequencies
disp('Calculated IF Frequencies for 6 NB-IoT Devices:')
disp(IF_carrier_freqs)
% Cell array to store waveform and sampling rate for each device
waveforms = cell(1, numDevices);
samplingRates = cell(1, numDevices);

% Loop through each device and generate waveform
for deviceIdx = 1:numDevices
    % Call the generate_NBiot_UL function
    [waveform, Fs] = Generate_NBiot();

    % Store waveform and sampling rate in cell arrays
    waveforms{deviceIdx} = waveform;
    samplingRates{deviceIdx} = Fs;
   
end
% Display information about the generated waveform for each device
disp(['Generated NB-IoT Uplink Waveform Information for Device ' num2str(deviceIdx) ':']);
disp(['Sampling Rate (Fs): ' num2str(Fs) ' Hz']);
disp(['Waveform Length: ' num2str(length(waveform)) ' samples']);
[signal,f] = pwelch(waveform,2048,1024,2048,Fs,"centered");
signal = signal/max(abs(signal));
generated = dsp.SpectrumAnalyzer('SampleRate', Fs, 'Title', 'The uplink generated waveform signal');
generated(waveform);
figure;plot(f/1e3,20*log10((signal)))
xlabel('Frequency in (KHz)')
ylabel('Normalized Magnitude (dB)')

% Create the composite IF signal
composite_IF_signal_upsampled = zeros(1,2*length(waveforms{1}));  % Initialize with the size of one waveform
IF_sampling_rate = 2 * max(IF_center_frequency + (numDevices - 1) * IF_frequency_spacing, Fs);
upsampling_factor = IF_sampling_rate/Fs;
upsampled_signal=zeros(1,2*length(waveforms{1}));

% Design the low-pass filter
% cutoff_frequency = IF_frequency_spacing / (2 * upsampling_factor);
% filter_order = 64;
% lpf = fir1(filter_order, cutoff_frequency / (Fs/2));
% Upsampling factor and filter specifications
 % Adjust as needed
 % Design and apply low-pass filter

for i = 1:numDevices 
        upsampled_signal = upsample((waveforms{i})', upsampling_factor);
         % Apply low-pass filter
%         lpf = lowpass('Fp,Fst,Ap,Ast', cutoff_frequency, cutoff_frequency * 1.5, 0.1, 80, Fs * upsampling_factor);
        cutoff_frequency = IF_frequency_spacing / (2 * upsampling_factor);

        filtered_signal= lowpass(upsampled_signal,cutoff_frequency,IF_sampling_rate);
%         filtered_signal = filter(lpf, 1, upsampled_signal);
        
        % Calculate the carrier frequency for the current device
%         carrier_frequency = IF_center_frequency + (deviceIndex - 1) * IF_frequency_spacing;
        
        % Perform digital upconversion by multiplying with a complex exponential
        t = (0:length(filtered_signal)-1)/ Fs / upsampling_factor  ;
        upconverted_signal = filtered_signal .* exp(1i*2 * pi * IF_carrier_freqs(i) * t);
        
        % Add the upconverted signal to the composite IF signal
        composite_IF_signal_upsampled = composite_IF_signal_upsampled + upconverted_signal;
end

% Plot the magnitude spectrum of the composite IF signal
[IF_signal_spectrum, IF_frequencies] = pwelch(composite_IF_signal_upsampled, 2048, 1024, 2048, Fs * upsampling_factor,"centered");
IF_signal_spectrum = IF_signal_spectrum / max(abs(IF_signal_spectrum));

% Plot the magnitude spectrum
figure;
plot(IF_frequencies/1e3, 10*log10((IF_signal_spectrum)));  % Using log scale for better visualization
xlabel('Frequency (kHz)')
ylabel('Power/Frequency (dB/Hz)')
title('Power Spectral Density Estimate using pwelch')
grid on;

IF_freq = 500e3;
composite_IF_signal_upsampled=composite_IF_signal_upsampled.*exp(1i*2 * pi * IF_freq * t);
%receiver 
received_waveforms_down_converged = zeros(6,length(composite_IF_signal_upsampled));
received_downsampled = zeros(6,length(waveform'));
receivedfiltered_signal=zeros(6,length(composite_IF_signal_upsampled'));
rx_signal_spectrum= zeros(6,length(signal));
f2=zeros(6,2048);
for i=1:6
    received_waveforms_down_converged(i,:)=composite_IF_signal_upsampled.* exp(1i*2 * pi *(IF_freq + IF_carrier_freqs(i)) * t);
%     [rx_signal_spectrum(i,:),f2(i,:)] = pwelch( received_waveforms_down_converged(i,:),2048,1024,2048,IF_sampling_rate,"centered");
%     rx_signal_spectrum(i,:) = rx_signal_spectrum(i,:)/max(abs(rx_signal_spectrum(i,:)));
    receivedfiltered_signal(i,:) = lowpass(received_waveforms_down_converged(i,:), 95000, IF_sampling_rate);
    downsampling_factor = upsampling_factor;
    received_downsampled(i,:) = downsample(receivedfiltered_signal(i,:), downsampling_factor);
    
    
end
% A=composite_IF_signal_upsampled*.exp()
figure;
subplot(4,1, 1);
pwelch(composite_IF_signal_upsampled, 2048, 1024, 2048, IF_sampling_rate, 'centered');
title('Spectrum after Digital Downconversion');
subplot(4,1, 2);
pwelch(received_waveforms_down_converged(1,:), 2048, 1024, 2048, IF_sampling_rate, 'centered');
title('Spectrum after Digital Downconversion');
subplot(4, 1, 3);
pwelch(receivedfiltered_signal(1,:), 2048, 1024, 2048, IF_sampling_rate, 'centered');
title('Spectrum after Digital Downconversion and FIR Filtering');
subplot(4, 1, 4);
pwelch(received_downsampled(1,:), 2048, 1024, 2048, IF_sampling_rate/downsampling_factor, 'centered');
title('Spectrum after Digital Downconversion, FIR Filtering, and Downsampling');

sgtitle('Spectral Analysis at Different Processing Stages');
% plot(f2()/1e3, 10*log10(fftshift(rx_signal_spectrum(2,:))));  % Using log scale for better visualization
% xlabel('Frequency (kHz)')
% ylabel('Power/Frequency (dB/Hz)')
% title('Power Spectral Density Estimate using pwelch')
% grid on;
