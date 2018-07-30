function threshold_record(x,DeltaT,threshold,fileID)

% initial the condition
if x(1) < threshold
    under_threshold = true;
else
    under_threshold = false;
end

% track for the trigger point
for i=1:length(x)
    if (x(i) > threshold) && under_threshold
        fprintf(fileID, 'Triggered on:\n%f ms, %d th data.\n', i*DeltaT, i);
        under_threshold = ~under_threshold;
    elseif (x(i) < threshold) && ~under_threshold
        fprintf(fileID, 'Triggered off:\n%f ms, %d th data.\n',i*DeltaT, i);
        under_threshold = ~under_threshold;
    end 
end