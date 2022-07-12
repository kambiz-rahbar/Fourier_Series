clc
clear
close all

signal_length_in_seconds = 2*pi;
number_of_signal_samples_in_second = 4000;
number_of_samples_in_signal = floor(signal_length_in_seconds * number_of_signal_samples_in_second);
delta_time =  signal_length_in_seconds/(number_of_samples_in_signal-1);
time = 0:delta_time:signal_length_in_seconds;

signal = 0*time;
signal(2*number_of_samples_in_signal/4 : 3*number_of_samples_in_signal/4) = 1;
signal = reshape(signal, 1, []);

number_of_terms = 10;%round(sqrt(length(signal)));
disp_results = 1;
[approximated_signal, Fourier_Series_terms, A0, A, B] = Fourier_Series(signal, signal_length_in_seconds, number_of_signal_samples_in_second, number_of_terms, disp_results);

