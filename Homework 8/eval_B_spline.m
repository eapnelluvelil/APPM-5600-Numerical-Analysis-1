%{
    Function: eval_B_spline
    Inputs:
        center - center of B-spline
        t - data point to evaluate B-spline at
        h - spacing between B-spline centers
    Outputs:
        val - value of B-spline evaluated at t 
%}

function val = eval_B_spline(center, t, h)
    x_i_m_2 = center - 2*h;
    x_i_m_1 = center - h;
    x_i_p_1 = center + h;
    x_i_p_2 = center + 2*h;

    B_spline_support = [x_i_m_2; x_i_m_1; center; x_i_p_1; x_i_p_2];
    
    % Determine which sub-interval t falls into
    bin = discretize(t, B_spline_support, 'IncludedEdge', 'left');

    if bin == 1
        val = (t - x_i_m_2)^(3);
    elseif bin == 2
        val = h^(3) + 3*h^(2)*(t - x_i_m_1) + 3*h*(t - x_i_m_1)^(2) - 3*(t - x_i_m_1)^(3);
    elseif bin == 3
        val = h^(3) + 3*h^(2)*(x_i_p_1 - t) + 3*h*(x_i_p_1 - t)^(2) - 3*(x_i_p_1 - t)^(3);
    elseif bin == 4
        val = (x_i_p_2 - t)^(3);
    else
        val = 0;
    end

    val = val/(h^(3));
end