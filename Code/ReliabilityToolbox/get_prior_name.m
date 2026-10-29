function dist_name = get_prior_name(pshape)
% Converts Dynare prior shape codes to distribution names
    switch pshape
        case 1
            dist_name = 'Beta';
        case 2
            dist_name = 'Gamma';
        case 3
            dist_name = 'Normal';
        case 4
            dist_name = 'Inv. Gamma';
        case 5
            dist_name = 'Uniform';
        case 6
            dist_name = 'Inv. Gamma 2';
        case 8
            dist_name = 'Weibull';
        otherwise
            dist_name = 'Unknown';
    end
end