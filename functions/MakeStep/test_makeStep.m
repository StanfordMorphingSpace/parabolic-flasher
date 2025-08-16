clear
import matlab.unittest.TestCase
import matlab.unittest.constraints.IsEqualTo
import matlab.unittest.constraints.AbsoluteTolerance

eps = 1e-15;
testCase = TestCase.forInteractiveUse;

for test_idx=1:2
    for test_type=["U", "F"]
    data = load("MakeStepTestJacobian_"+test_type+"_"+num2str(test_idx)+".mat");
    p = data.p; 
    v = data.v;
    a = data.a;
    varargin = data.varargin{1};
    expected_p = data.p_new;
    expected_v = data.v_new;
    expected_E_v = data.E_v;
    labels = data.labels;
    beta = data.beta;
    i = data.i;
    mass = data.mass;
    dt = data.dt;
    getb = data.getb;
    if test_type == "U"
        getJacobian = @getJacobian_u;
        getJacobianFast = @getJacobian_u_fast;
    else
        getJacobian = @getJacobian_f;
        getJacobianFast = @getJacobian_f_fast;
    end

    [actual_p, actual_v, actual_E_v] = makeStep(p, v, a, getJacobian, getb, labels, beta, i, mass, dt, varargin);
    testCase.verifyThat(actual_p,IsEqualTo(expected_p, "Within",AbsoluteTolerance(eps)))
    testCase.verifyThat(actual_v,IsEqualTo(expected_v, "Within",AbsoluteTolerance(eps)))
    testCase.verifyThat(actual_E_v,IsEqualTo(expected_E_v, "Within",AbsoluteTolerance(eps)))


    [actual_p, actual_v, actual_E_v] = makeStep(p, v, a, getJacobianFast, getb, labels, beta, i, mass, dt, varargin);
    testCase.verifyThat(actual_p,IsEqualTo(expected_p, "Within",AbsoluteTolerance(eps)))
    testCase.verifyThat(actual_v,IsEqualTo(expected_v, "Within",AbsoluteTolerance(eps)))
    testCase.verifyThat(actual_E_v,IsEqualTo(expected_E_v, "Within",AbsoluteTolerance(eps)))
    end
end