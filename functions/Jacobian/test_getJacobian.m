clear
import matlab.unittest.TestCase
import matlab.unittest.constraints.IsEqualTo
import matlab.unittest.constraints.AbsoluteTolerance
eps=1e-8;
testCase = TestCase.forInteractiveUse;

for i=1:2
    data = load("TestGetJacobian_U_"+num2str(i)+".mat");
    beta = data.beta; 
    labels = data.labels;
    nodes_f = data.nodes_f;
    varargin = data.varargin;
    expectedJ = data.J;
    
    actualJ = getJacobian_u(nodes_f, labels, beta, varargin);
    % assert(isequal(actualJ, expectedJ), "Test "+num2str(i)+" failed: Jacobian does not match expected output.")
    testCase.verifyThat(actualJ,IsEqualTo(expectedJ, "Within",AbsoluteTolerance(eps)))

    actualJ = getJacobian_u_fast(nodes_f, labels, beta, varargin);
    % assert(isequal(actualJ, expectedJ), "Test "+num2str(i)+" failed: Jacobian fast does not match expected output.")
    testCase.verifyThat(actualJ,IsEqualTo(sparse(expectedJ), "Within",AbsoluteTolerance(eps)))

end