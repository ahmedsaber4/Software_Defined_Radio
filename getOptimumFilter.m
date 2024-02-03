

function optimum_filter = getOptimumFilter(filters, orders, passbandRipples, stopbandAttenuations, transitionBW)
    scores = zeros(1, length(filters));    
    for i = 1 : length(filters)
        scores(i) = evaluateFilter(orders(i), passbandRipples(i), stopbandAttenuations(i), transitionBW(i));
    end
    [~, idx] = max(scores);
    disp(scores);
    optimum_filter = filters(idx);
end