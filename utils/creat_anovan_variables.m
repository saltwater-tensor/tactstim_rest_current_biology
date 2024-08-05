function [Y,Y_state,Y_response] = creat_anovan_variables(temporal_metric,response_type)
Y = [];
Y_state = [];
Y_response = [];

for st = 1:size(temporal_metric,2)
    state_variable = ['State' num2str(st)];

    for trials = 1:size(temporal_metric,1)

        try
            Y = [Y;temporal_metric{trials,st}];

            ystate = cell((length(temporal_metric{trials,st})),1);
            ystate(:) = {state_variable};
            y_resp = cell((length(temporal_metric{trials,st})),1);
            y_resp(:) = {response_type};

            Y_state = [Y_state;ystate];
            Y_response = [Y_response;y_resp];

        catch

            Y = [Y;temporal_metric(trials,st)];

            ystate = cell((length(temporal_metric(trials,st))),1);
            ystate(:) = {state_variable};
            y_resp = cell((length(temporal_metric(trials,st))),1);
            y_resp(:) = {response_type};

            Y_state = [Y_state;ystate];
            Y_response = [Y_response;y_resp];

        end



    end
end





end