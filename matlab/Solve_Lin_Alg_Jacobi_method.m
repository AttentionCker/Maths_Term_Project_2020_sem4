fprintf('\nTake Ax = b as the system of linear equations in x\n');
A = input("Enter The COEFFICIENT Matrix: A = ");
b = input("Enter The COLUMN vector b = ");


if ( size(A,1) ~= size(A,2) ) 
    error('Entered Matrix A is not a square Matrix');

else    % finding the solution upto the desired tolerance 
    
    tol = input("Enter the desired tolerance for the maximum acceptable error: ");
    
    % Setting initial values for solving
    n = size(A,1);

    % Assuming an intial value
    x0 = zeros(n,1);

    % Initializing values for D and L+U
    D = diag(diag(A));
    L_plus_U = A - D;  

    % epsilon set for max error
    epsilon = [Inf; zeros(n-1, 1)];
    maxIterations = 500;
    itr = 1;
    while ( max(epsilon) > tol && itr < maxIterations )

        % Computing x using Jacobi method
        x = D\(b - L_plus_U*x0);        % (matrix left division) D\N is equivalent to multiplying by inv(D) ... where D is a diagonal matrix  
        
        % Computing error vector epsilon
        epsilon = abs(x - x0);

        % updating
        x0 = x;
        itr = itr + 1;
    end

    fprintf('\n\nThe solution computed is \nx = \n');
    for i = 1: n
        fprintf("%.10f\n", x(i));
    end    
end