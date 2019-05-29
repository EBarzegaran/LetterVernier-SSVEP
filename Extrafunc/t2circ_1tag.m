function [Z_est, confidence_radii, p, t2circ] = t2circ_1tag(ft_data,alpha)
    % T^2_circ statistic based on Victor & Mast (1991)
    % Usage: [Z_est critical_range_zero_center p t2circ] = t2circ_1tag(ft_data,alpha)
    % Inputs: 
    % ft_data- Complex[M x N]: M electrodes by N trials of complex valued fourier coefficients
    % alpha- significance value for determining confidence interval: ie. 0.05
    % Outputs:
    % Z_est- Mean fourier
    % confidence_radii- N-length array of real numbers representing the 2D confidence circle radius  for each frequency
    % p- M-length array of significance values for each frequency
    % t2circ- M-length array of t2circ for each frequency

    % check inputs
    if nargin < 1
        disp('Using defaults: random data 10sweeps of 1000 samples')
        ft_data = fft(rand(10,1000));
    end
    if nargin < 2
        alpha = 0.05;
    else
    end

    if alpha<0 || alpha>1 
        error('alpha must be between 0 and 1');
    end

    if isreal(ft_data)
        error('ft_data are not complex: are the FT-ed?');
    end

    % determine number of sweeps and number of frequencies
    M = size(ft_data,1);
    N = size(ft_data,2);

    % compute estimate of population mean <Z>est from V & M 1991
    Z_est = mean(ft_data,2);

    % compute V_indiv (Equation 1: V & M 1991)
    %First compute the squared differences from the mean for each trial- Variance
    freq_diff = abs(ft_data-repmat(Z_est,1,N)).^2;

    % Sum the squared differences
    freq_diff_sum = sum(freq_diff,2);

    % Normalize by number of trials
    v_indiv = 1/(2*(N -1)) .* freq_diff_sum;

    % Compute V_Group (Equation 2: V & M 1991)
    % NOTE for almost all purposes population mean is assumed to be zero as it
    % is here. If for whatever reason you want to test agains something
    % non-zero you will need to modify the code to include an N-length array of
    % complex valued poluation means and modify v_group beloq

    v_group =  (abs(Z_est).^2) * N/2;

    % compute the tcirc-statistic
    t2circ = (v_group./v_indiv)/N;

    % perform an F-test/Compute Confidence interval
    % Find critical F for corresponding alpha level drawn from F-distribution F(2,2N-2)
    f_crit = finv(1-alpha,2,2*N-2);

    % NOTE: it is NxT2Circ that is distributed according to F(2,2N-2)
    p = fpdf(N*t2circ,2,2*N- 2);

    %Set bounds of confidence interval [equations 4/5]
    confidence_radii = sqrt((f_crit/N) * 2*v_indiv);
end




