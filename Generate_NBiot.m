

function [waveform,Fs]= Generate_NBiot()
tbs = 144;                          % The transport block size
totNumBlks = 10;                     % Number of simulated transport blocks
subspacing='15kHz';
ue = struct();                      % Initialize the UE structure
ue.NBULSubcarrierSpacing = subspacing ; % 3.75kHz, 15kHz
ue.NNCellID = 0;                    % Narrowband cell identity

chs = struct();
% NPUSCH carries data or control information
chs.NPUSCHFormat = 'Data'; % Payload type (Data or Control)
% The number of subcarriers used for NPUSCH 'NscRU' depends on the NPUSCH
% format and subcarrier spacing 'NBULSubcarrierSpacing' as shown in TS
% 36.211 Table 10.1.2.3-1. There are 1,3,6 or 12 contiguous subcarriers for
% NPUSCH
chs.NBULSubcarrierSet = 0:11;  % Range is 0-11 (15kHz); 0-47 (3.75kHz)
%chs.NRUsc = length(chs.NBULSubcarrierSet);
chs.NRUsc = 12;
chs.CyclicShift = 0;   % Cyclic shift required when NRUsc = 3 or 6
chs.RNTI = 0;          % RNTI value
chs.NLayers = 1;       % Number of layers
chs.NRU = 2;           % Number of resource units
chs.NRep = 4;          % Number of repetitions of the NPUSCH
chs.SlotIdx = 0;       % Start slot index in a bundle
% The symbol modulation depends on the NPUSCH format and NscRU as
% given by TS 36.211 Table 10.1.3.2-1
chs.Modulation = 'QPSK';
rvDCI = 0;             % RV offset signaled via DCI (See 36.213 16.5.1.2)

% Specify the NPUSCH and DM-RS power scaling in dB for plot visualization
chs.NPUSCHPower = 30;
chs.NPUSCHDRSPower = 34;

chs.SeqGroupHopping = 'on'; % Enable/Disable Sequence-Group Hopping for UE
chs.SeqGroup = 0;           % Delta_SS. Higher-layer parameter groupAssignmentNPUSCH

% Get number of time slots in a resource unit NULSlots according to
% TS 36.211 Table 10.1.2.3-1
if strcmpi(chs.NPUSCHFormat,'Data')
    if chs.NRUsc == 1
        %NULSlots = 16;
        NULSlots = 16;
    elseif any(chs.NRUsc == [3 6 12])
       NULSlots = 24/chs.NRUsc;
       % NULSlots = 2;

    else
        error('Invalid number of subcarriers. NRUsc must be one of 1,3,6,12');
    end
elseif strcmpi(chs.NPUSCHFormat,'Control')
    %NULSlots = 4;
            NULSlots = 4;

else
    error('Invalid NPUSCH Format (%s). NPUSCHFormat must be ''Data'' or ''Control''',chs.NPUSCHFormat);
end
chs.NULSlots = NULSlots;

NSlotsPerBundle = chs.NRU*chs.NULSlots*chs.NRep; % Number of slots in a codeword bundle
TotNSlots = totNumBlks*NSlotsPerBundle;   % Total number of simulated slots



% Initialize the random generator to default state
rng('default');

% Get the slot grid and number of slots per frame
emptySlotGrid = lteNBResourceGrid(ue);
slotGridSize = size(emptySlotGrid);
NSlotsPerFrame = 20/(slotGridSize(1)/12);

state = [];    % NPUSCH encoder and DM-RS state, auto re-initialization in the function
trblk = [];    % Initialize the transport block
txgrid = [];   % Full grid initialization

% Display the number of slots being generated
% fprintf('\nGenerating %d slots corresponding to %d transport block(s)\n',TotNSlots,totNumBlks);
for slotIdx = 0+(0:TotNSlots-1)
    % Calculate the frame number and slot number within the frame
    ue.NFrame = fix(slotIdx/NSlotsPerFrame);
    ue.NSlot = mod(slotIdx,NSlotsPerFrame);

    if isempty(trblk)
       if strcmpi(chs.NPUSCHFormat,'Data')
           % UL-SCH encoding is done for the two RV values used for
           % transmitting the codewords. The RV sequence used is determined
           % from the rvDCI value signaled in the DCI and alternates
           % between 0 and 2 as given in TS 36.213 Section 16.5.1.2

           % Define the transport block which will be encoded to create the
           % codewords for different RV
           trblk = randi([0 1],tbs,1);

           % Determine the coded transport block size
           [~, info] = lteNPUSCHIndices(ue,chs);
           outblklen = info.G;
           % Create the codewords corresponding to the two RV values used
           % in the first and second block, this will be repeated till all
           % blocks are transmitted
           chs.RV = 2*mod(rvDCI+0,2); % RV for the first block
           cw = lteNULSCH(chs,outblklen,trblk); % CRC and Turbo coding
           chs.RV = 2*mod(rvDCI+1,2); % RV for the second block
           cw = [cw lteNULSCH(chs,outblklen,trblk)]; %#ok<AGROW> % CRC and Turbo coding is repeated
       else
           trblk = randi([0 1],1); % 1 bit ACK
           % For ACK, the same codeword is transmitted every block as
           % defined in TS 36.212 Section 6.3.3
           cw = lteNULSCH(trblk);
       end
       blockIdx = 0; % First block to be transmitted
    end

    % Initialize grid
    slotGrid = emptySlotGrid;

    % NPUSCH encoding and mapping onto the slot grid
    txsym = lteNPUSCH(ue,chs,cw(:,mod(blockIdx,size(cw,2))+1),state);

    % Map NPUSCH symbols in the grid of a slot
    indicesNPUSCH = lteNPUSCHIndices(ue,chs);
    slotGrid(indicesNPUSCH) = txsym*db2mag(chs.NPUSCHPower);

    % Create DM-RS sequence and map to the slot grid
    [dmrs,state] = lteNPUSCHDRS(ue,chs,state);
    indicesDMRS = lteNPUSCHDRSIndices(ue,chs);
    slotGrid(indicesDMRS) = dmrs*db2mag(chs.NPUSCHDRSPower);

    % Concatenate this slot to the slot grid
    txgrid = [txgrid slotGrid]; %#ok<AGROW>

    % If a full block is transmitted, increment the clock counter so that
    % the correct codeword can be selected
    if state.EndOfBlk
        blockIdx = blockIdx + 1;
    end

    % trblk err count and re-initialization
    if state.EndOfTx
       % Re-initialize to enable the transmission of a new transport block
       trblk = [];
    end

end

% Perform SC-FDMA modulation to create time domain waveform
ue.CyclicPrefixUL = 'Normal'; % Normal cyclic prefix length for NB-IoT
[waveform,scfdmaInfo] = lteSCFDMAModulate(ue,chs,txgrid);

% Create an image of overall resource grid
% figure
% im = image(abs(txgrid));
% cmap = parula(64);
% colormap(im.Parent,cmap);
% axis xy;
% title(sprintf('NB-IoT Uplink RE Grid (NRep = %d, NRUsc = %d, NRU = %d)',chs.NRep,chs.NRUsc,chs.NRU))
% xlabel('OFDM symbols')
% ylabel('Subcarriers')
% % Create the legend box to indicate the channel/signal types associated with the REs
% reNames = {'NPUSCH';'DM-RS'};
% clevels = round(db2mag([chs.NPUSCHPower chs.NPUSCHDRSPower]));
% N = numel(reNames);
% L = line(ones(N),ones(N), 'LineWidth',8); % Generate lines
% % Set the colors according to cmap
% set(L,{'color'},mat2cell(cmap( min(1+clevels,length(cmap) ),:),ones(1,N),3));
% legend(reNames{:});
Fs= scfdmaInfo.SamplingRate;
end