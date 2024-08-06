function outputString = getPathSubset(path,nFromEnd)
    splitPath = split('/',path);
    outputString = [];
    for i = 1:length(splitPath) - nFromEnd
        outputString = [outputString, '/', splitPath{i}];
    end
end