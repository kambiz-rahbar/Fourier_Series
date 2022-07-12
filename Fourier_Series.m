function [approximated_signal, Fourier_Series_terms, A0, A, B] = Fourier_Series(signal, signal_length_in_seconds, number_of_signal_samples_in_second, number_of_terms, disp_results)
% https://en.wikipedia.org/wiki/Fourier_series
% https://www.youtube.com/watch?v=PP5ox7evg7o&t=72s

number_of_samples_in_signal = floor(signal_length_in_seconds * number_of_signal_samples_in_second);
delta_time =  signal_length_in_seconds/(number_of_samples_in_signal-1);
time = 0:delta_time:signal_length_in_seconds;

A0 = sum(signal.*ones(size(time)))*delta_time/pi;

Fourier_Series_terms = zeros(number_of_terms, size(signal,2));
Fourier_Series_terms(1, :) = A0/2*ones(size(signal));
A = zeros(number_of_terms-1, 1);
B = zeros(number_of_terms-1, 1);
for k = 1:number_of_terms - 1
    A(k) = sum(signal.*cos(2*pi*k*time/signal_length_in_seconds))*delta_time/pi;
    B(k) = sum(signal.*sin(2*pi*k*time/signal_length_in_seconds))*delta_time/pi;

    Fourier_Series_terms(k+1, :) = A(k)*cos(2*k*pi*time/signal_length_in_seconds) + B(k)*sin(2*k*pi*time/signal_length_in_seconds);
end
approximated_signal = sum(Fourier_Series_terms, 1);

if disp_results
    figure(disp_results);

    subplot(2, 2 ,1);
    plot(time, Fourier_Series_terms);
    grid minor;
    xlabel('time');
    ylabel('magnitude');
    title('terms of fourier series');
    xlim([min(time),max(time)]);

    subplot(2, 2 ,2);
    plot(time, signal), hold;
    plot(time, approximated_signal);
    grid minor;
    xlabel('time');
    ylabel('magnitude');
    legend('signal','aproximated signal');
    title('fourier series aproximation');
    xlim([min(time),max(time)]);

    acumulated_sum = 0;
    aproximation_error = zeros(number_of_terms, 1);
    for k = 1:number_of_terms
        acumulated_sum = acumulated_sum + Fourier_Series_terms(k, :);
        aproximation_error(k) = norm(signal-acumulated_sum)./norm(signal);
    end
    subplot(2, 2 ,3);
    semilogy(aproximation_error);
    grid minor;
    xlabel('up to term');
    ylabel('error');
    title('aproximation error');
    xlim([0, number_of_terms]);
    
end