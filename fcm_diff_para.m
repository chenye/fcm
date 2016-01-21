
function [avg] = fcm_diff_para(filename)

avg = zeros(9, 9, 7);

result = importdata(filename);
for i=1:9
    for j=1:9
        line_id_start = ((i-1)*9+j-1)*20+1;
        line_id_end = line_id_start + 19;
        if line_id_end > size(result,1)
            break;
        end
        avg(i,j,:) = mean(result( line_id_start : line_id_end, : ));
    end
end
