%Returns the individual compatibility squared chi2 distance for the hypothesis
%Measurement is 2 x 1, range bearing measurement.
%Hypothesis is 1 x 1 representing index of expected landmark.
function result = individual_compatibility(measurement,hypothesis)
    result = joint_compatibility(measurement,hypothesis);
end