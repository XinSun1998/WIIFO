%% Field II: 16-by-16 2D array with user-defined delays and apodization
clear;
clc;
close all;

%% ============================================================
% 1. Field II initialization
addpath('FieldII Path')
% =============================================================
field_init(-1);

Lx = -1; % 1 2.5 1.5
Rx = 1;  % 1 2.5 1.5
Ly = -1; % 1 2.5 1.5
Ry = 1;  % 1 2.5 1.5


c  = 1480;             % Speed of sound [m/s]
f0 = 10e6;             % Center frequency [Hz]
fs = 200e6;            % Field II sampling frequency [Hz]

set_field('c', c);
set_sampling(fs);

lambda = c / f0;

addpath('YourLocalPath')

Filename = 'SPIRAL_';
load([Filename,'pattern_7mm_correctedV.mat'], 'elementDelayNsMatrix', ...
                      'elementAmplitudeMatrix','elementPhaseMatrix');

load('Patterncmap_5.mat')
load('Phasecmap_light.mat')
load('Apodcmap.mat')
%% ============================================================
% 2. Array parameters
% =============================================================
Nx = 16;
Ny = 16;

pitch_x = 300e-6;      % Element center-to-center pitch [m]
pitch_y = 300e-6;

kerf_x = 20e-6;
kerf_y = 20e-6;

element_width  = pitch_x - kerf_x;   % 270 um
element_height = pitch_y - kerf_y;   % 270 um

sub_x = 1;
sub_y = 1;

% initial_focus = [0 0 10e-3];
initial_focus = [0 0 0];
% Enable all physical elements
enabled = ones(Ny, Nx);

Th = xdc_2d_array( ...
    Nx, ...
    Ny, ...
    element_width, ...
    element_height, ...
    kerf_x, ...
    kerf_y, ...
    enabled, ... 
    sub_x, ...
    sub_y, ...
    initial_focus);

%% ============================================================
% 3. Excitation waveform and element impulse response
% =============================================================

nCycles = 1;
tExc = 0 : 1/fs : nCycles/f0;
excitation = sin(2*pi*f0*tExc);

% Band-limited element impulse response
fractionalBandwidth = 0.60;

tImp = -2/f0 : 1/fs : 2/f0;

impulseResponse = gauspuls( ...
    tImp, ...
    f0, ...
    fractionalBandwidth);

xdc_excitation(Th, excitation);
xdc_impulse(Th, impulseResponse);

%% ============================================================
% 4. Input  delay and apodization matrices here
% =============================================================

delayMatrixNs = elementDelayNsMatrix;
apodMatrix    = elementAmplitudeMatrix;
% -------------------------------------------------------------

assert(isequal(size(delayMatrixNs), [Ny, Nx]), ...
    'delayMatrixNs must be a 16-by-16 matrix.');

assert(isequal(size(apodMatrix), [Ny, Nx]), ...
    'apodMatrix must be a 16-by-16 matrix.');

assert(all(isfinite(delayMatrixNs(:))), ...
    'delayMatrixNs contains NaN or Inf.');

assert(all(isfinite(apodMatrix(:))), ...
    'apodMatrix contains NaN or Inf.');

%% ============================================================
% 5. Convert matrices to Field II element vectors
% =============================================================

% Field II internally uses a one-dimensional vector for the
% physical elements.
%
% For the array definition used here, start with MATLAB column-major
% ordering. The orientation can be checked using the single-element
% test included near the end of this script.

delayVector = reshape(delayMatrixNs, 1, []);
apodVector  = reshape(apodMatrix,    1, []);

% Convert nanoseconds to seconds
delayVectorSec = delayVector * 1e-9;

% Shift the delays so the minimum delay is zero.
% This does not change the relative phase pattern.
delayVectorSec = delayVectorSec - min(delayVectorSec);

% Normalize apodization without changing its relative values
maxApod = max(abs(apodVector));

if maxApod > 0
    apodVector = apodVector / maxApod;
else
    error('All apodization coefficients are zero.');
end

% Apply user-defined delays
xdc_focus_times(Th, 0, delayVectorSec);

% Apply user-defined apodization
xdc_apodization(Th, 0, apodVector);

%% ============================================================
% 6. Display delay and apodization matrices
% =============================================================
xElement = ((0:Nx-1) - (Nx-1)/2) * pitch_x * 1e3;
yElement = ((0:Ny-1) - (Ny-1)/2) * pitch_y * 1e3;

figure
imagesc(xElement, yElement, apodMatrix);
axis image;

% ===== Remove ALL axes/ticks =====
ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.LineWidth = 2;   

set(gca, 'YDir', 'normal');
title('ElementApodMatrix')
colormap(Apodcmap)


figure
imagesc(xElement, yElement, elementPhaseMatrix);
axis image;

% ===== Remove ALL axes/ticks =====
ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.LineWidth = 2;   

set(gca, 'YDir', 'normal');
title('ElementPhaseMatrix')
% colorbar
clim([-pi pi]);
colormap(Phasecmap)


%% ============================================================
% 7. Calculate X-Y pressure field at a selected depth
% =============================================================

observationDepth = 7e-3;

xXY = linspace(-3e-3, 3e-3, 301);
yXY = linspace(-3e-3, 3e-3, 301);

pressureXY = zeros(length(yXY), length(xXY));

fprintf('\nCalculating X-Y field at z = %.2f mm...\n', ...
    observationDepth * 1e3);

for iy = 1:length(yXY)

    fprintf('X-Y plane: %d / %d\n', iy, length(yXY));

    for ix = 1:length(xXY)

        fieldPoint = [xXY(ix), yXY(iy), observationDepth];

        [hp, ~] = calc_hp(Th, fieldPoint);

        if isempty(hp)
            pressureXY(iy, ix) = 0;
        else
            pressureXY(iy, ix) = sqrt(sum(abs(hp(:)).^2));
        end
    end
end

pressureXY = pressureXY / max(pressureXY(:));

pressureXYdB = 20 * log10(pressureXY + eps);
pressureXYdB(pressureXYdB < -50) = -50;

%% Plot X-Y field

figure('Name', 'X-Y pressure field');

imagesc( ...
    xXY * 1e3, ...
    yXY * 1e3, ...
    pressureXYdB');

axis image;
set(gca, 'YDir', 'normal');
ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.LineWidth = 2;      
ax.Box = 'on';          

clim([-20 0]);
xlim([Lx Rx])
ylim([Ly,Ry])
colormap(Patterncmap)
title('SimulatedPattern')

%% ============================================================
% 8. Free Field II aperture
% =============================================================

xdc_free(Th);
field_end;

fprintf('\nSimulation completed.\n');