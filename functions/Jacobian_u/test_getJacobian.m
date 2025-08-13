clear
for i=1:2
    data = load("Test_"+num2str(i)+".mat");
    beta = data.beta; 
    labels = data.labels;
    nodes_f = data.nodes_f;
    varargin = data.varargin{1};
    expectedJ = data.J;
    
    actualJ = getJacobian_u(nodes_f, labels, beta, varargin);
    assert(isequal(actualJ, expectedJ), "Test "+num2str(i)+" failed: Jacobian does not match expected output.")

    actualJ = getJacobian_u_fast(nodes_f, labels, beta, varargin);
    assert(isequal(actualJ, expectedJ), "Test "+num2str(i)+" failed: Jacobian fast does not match expected output.")
end
disp("All tests passed successfully.");