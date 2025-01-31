% Define the template file name and target file names
templateFile = 'sAIF_BirthDeath_mu1.m';

for i = 2:32
    % Read the content of the template file
    fileContent = fileread(templateFile);
    
    % Replace Script_Index = 1 with the desired value
    newContent = regexprep(fileContent, 'Script_Index = 1;', sprintf('Script_Index = %d;', i));
    
    % Define the new file name
    newFileName = sprintf('sAIF_BirthDeath_mu%d.m', i);
    
    % Write the new content to the new file
    fid = fopen(newFileName, 'w');
    fwrite(fid, newContent);
    fclose(fid);
end
