function total_score = evaluateFilter(order, passbandRipple, stopbandAttenuation, transitionBandwidth)    
    param = [order, passbandRipple, stopbandAttenuation, transitionBandwidth];
    % Evaluation criteria and corresponding scores
    orderCriteria = [300, 200, 100];
    rippleCriteria = [0.4, 0.3, 0.2];
    attenuationCriteria = [1/-15, 1/-25, 1/-60];
    bandwidthCriteria = [80, 50, 30];
    
    criteria = [orderCriteria; rippleCriteria; attenuationCriteria; bandwidthCriteria];
    total_score = 0;

    for i = 1 : length(criteria)
        score = find(param(i) >= criteria(i,:), 1, 'first');
        if isempty(score)
            score = length(criteria(i,:));
        else
            score = score-1;
        end
        total_score = total_score + score;
    end
end