function R_hat = calculate_reconstruction_metric(Datastream,Dictionary,A,lambda)
%calculate_reconstruction_metric uses the following formula to assign a
%score to the reconstruction:
%   (1/2)||x - D*a||_2^2 + lambda*||a||_1
%or, extended to all datapoints,
%   (1/2)||X - D*A||_2^2 + lambda*||A||_1

Reconstruction = Dictionary * A;
term_1 = 0.5 * vecnorm(Datastream - Reconstruction, 2, 1).^2;
term_2 = lambda * vecnorm(A, 1, 1);

R_hat = term_1 + term_2;

end

