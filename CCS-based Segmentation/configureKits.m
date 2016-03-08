function statString = configureKits(statNum, msg)
% 用途：
% 1）进行轮廓演化前，使用该函数进行环境配置，关闭所有窗口，清除命令行内容，然后添加各个Toolkit的路径；
% 2）在轮廓演化算法调试时，通过调用本函数，实现对轮廓演化工具包的实时版本控制与远程代码提交
% 要求：用户如需使用本函数提供的版本控制功能，需要首先在本地下载Git的客户端（链接：http://git-scm.com/download/win）
% 下载完成后，确认Git可执行文件的路径已经被添加到了系统路径（Path）中，一般而言安装时是会默认添加的
% 此后在正确配置了Git的基础上，才能调用本函数
% 
% 输入参数
% statNum：int类型变量，函数调用码
% 当statNum为0，表示调用该函数完成轮廓演化的基本环境配置
% 当statNum为1，表示调用该函数从Git的远程代码仓库（即Git服务器端）
% 当statNum为2，表示调用该函数将已有的修改全部添加到Git的本地暂存区
% 当statNum为3，表示调用该函数将已有的修改全部提交到Git的本地代码仓库
% 当statNum为4，表示调用该函数将所有的修改全部提交到Git的远程代码仓库
% msg：字符串类型变量，当statNum为3时，最好使用msg输入对当前所有修改的注释；如用户没有输入msg，则函数会提交默认注释
% 输出参数
% statString：字符串类型变量，configureKits函数的执行结果
%
% 局限：函数目前仅适用于Windows操作系统，Linux版本尚未开发

%% 基本参数处理
try
    statusNum = statNum;
catch  % 如果用户没有输入状态码
    statusNum = 0; % 默认为0
end

%% 当statusNum为0，该函数用于配置环境
if statusNum == 0
    close all
    clc;
    format long g
    addpath('./CCSToolkit'); % 添加CCS Toolkit工具包
    addpath('./EvolveToolkit'); % 添加Evolve Toolkit工具包
    addpath('./CCSDebugKit'); % 添加CCS演化的Debug工具包
    statString = 'Environment configuration has been completed';

%% 当statusNum为1，从Git的远程仓库下载代码到本地用于更新
elseif statusNum == 1
    [statCode, resultString] = system('git pull');
    if statCode == 0
         % 如果返回的函数执行结果字符串中包含'Already up-to-date',说明本地代码的版本已经是Git服务器端的最新版本
         % 或者比Git服务器端的版本还要新，也即存在本地修改还没提交到服务器端的代码
        if ~isempty(strfind(resultString, 'Already up-to-date'))
            statString = 'The code in this computer is already the newest version';
        % 如果返回的函数执行结果字符串中不包含'Already up-to-date',则说明本地需要更新，且已通过git
        % pull命令更新到了最新版本
        else
            statString = 'The code in this computer has been updated to the newest version';
        end
    else
        if ~isempty(strfind(resultString, 'Please, commit your changes or stash them before you can merge.'))
            statString = 'You should use configureKits(2) and configureKits(3) to commit your changes to the local repository before pulling';
        elseif ~isempty(strfind(resultString, 'Automatic merge failed; fix conflicts and then commit the result.'))
            statString = 'Successfully pulling from remote repository of Git but conflicts occurs, please check the conflicts and fix them';
        else
            statString = 'Failed to get the newest version of code from remote repository of Git';
        end       
    end

%% 当statusNum为2，提交修改到Git的本地暂存区    
elseif statusNum == 2
    [~, resultString] = system('git status');
    checkStr1 = 'Changes not staged for commit:';
    checkStr2 = 'Changes to be committed:';
    checkStr3 = 'nothing to commit, working directory clean';
    checkStr4 = 'You have unmerged paths.';
    checkStr5 = 'nothing added to commit but untracked files present';
    % 如果返回的函数执行结果字符串中包含checkStr1的内容,说明有修改尚未提交到暂存区，此时才能调用git add -A
    if ~isempty(strfind(resultString, checkStr1)) 
        [statCode, ~] = system('git add -A');
        if statCode == 0 % 如果状态码为0，表明正确提交
            statString = 'The changes of code have been added to the local stage area';
        else
            statString = 'Failed to add the changes to the local stage area';
        end
    elseif ~isempty(strfind(resultString, checkStr4))
        statString = 'Perhaps You have fixed all the conflicts but are still not able to commit, try to make some modification and then use "configureKits(2)" again';
    elseif ~isempty(strfind(resultString, checkStr2)) 
        statString = 'Nothing needed to be added to the local stage area. Please use "configureKits(3)" to commit changes to the local repository of Git ';
    elseif ~isempty(strfind(resultString, checkStr3))
        statString = 'Nothing needed to be committed';
    elseif ~isempty(strfind(resultString, checkStr5))
        [statCode, ~] = system('git add -A');
        if statCode == 0 % 如果状态码为0，表明正确提交
            statString = 'New Files detected and added to the local stage area';
        else
            statString = 'Failed to add the new files to the local stage area';
        end
    else
        statString = 'Failed to add the changes to the local stage area for being unable to detect the Git status'; 
    end
    
%% 当statusNum为3，提交修改到Git的本地仓库
elseif statusNum == 3
    try
        msgString = msg;
    catch  % 如果用户没有输入将代码提交到本地仓库时的注释
        msgString = 'new modifications'; % 默认注释为‘new modifications’
    end
    cmdString = ['git commit -m "' msgString '"'];
    [statCode, resultString] = system(cmdString);
    if statCode == 0 % 如果状态码为0，表明正确提交
        statString = 'The changes of code have been committed to the local repository of Git';
    else % 如果状态码不为0
        checkStr1 = 'nothing to commit, working directory clean';
        checkStr2 = 'Changes not staged for commit:';
        checkStr3 = 'fatal: Exiting because of an unresolved conflict.';
        % 如果返回的函数执行结果字符串中包含checkStr1的内容,说明是无新内容提交，而非提交错误
        if ~isempty(strfind(resultString, checkStr1)) 
            statString = 'Nothing new to commit to the local repository of Git, working directory is clean';
        % 如果返回的函数执行结果字符串中包含checkStr2的内容,说明用户尚未调用git add命令提交修改到暂存区，此时需要先调用configureKits(2)
        elseif ~isempty(strfind(resultString, checkStr2)) 
            statString = 'Failed to commit the changes, please use "configureKits(2)" to add the changes to the stage area first';
        elseif ~isempty(strfind(resultString, checkStr3)) 
            statString = 'Failed to commit the changes, please use "configureKits(2)" to fix the conflicts of merge first';
        else % 如果不包含checkStr1和checkStr2的内容，说明提交错误
            statString = 'Failed to commit the changes to the local repository of Git for unknown reason';
        end
    end

%% 当statusNum为4，提交修改到Git的远程仓库
elseif statusNum == 4
    % 由于Git自身不会检查用户在push之前是否已经在本地add和commit，因此在提交代码到远程仓库前，需要首先检查此前的修改是否都已经提交
    [~, resultString] = system('git status');
    checkStr1 = 'Changes not staged for commit:';
    checkStr2 = 'Changes to be committed:';
    checkStr3 = 'nothing to commit, working directory clean';
    % 如果返回的函数执行结果字符串中包含checkStr1的内容,说明还有若干修改尚未添加到暂存区
    if ~isempty(strfind(resultString, checkStr1))
        statString = 'Failed to commit the changes, please use "configureKits(2)" to add the changes to the stage area first';
    % 如果返回的函数执行结果字符串中包含checkStr2的内容,说明所有的修改均已添加到暂存区，但还尚未提交到本地仓库
    elseif ~isempty(strfind(resultString, checkStr2))
        statString = 'Failed to commit the changes, please use "configureKits(3, your notes)" to add the changes to the local repository first';
    % 如果返回的函数执行结果字符串中包含checkStr3的内容,说明所有的修改已提交到本地仓库
    elseif ~isempty(strfind(resultString, checkStr3))
    % 通过了上述验证，才会调用git pull进行代码的远程提交
        [statCode, resultString] = system('git push');
        checkStr4 = 'Everything up-to-date';
        if statCode == 0
            if ~isempty(strfind(resultString, checkStr4))
                statString = 'No need for pushing, for the code of the remote repository of Git is already the newest version';
            else
                statString = 'The changes of code have been committed to the remote repository of Git';
            end
        else
            statString = 'Failed to commit the changes to the remote repository of Git for unknown reason'; 
        end
    else
        statString = 'Failed to commit the changes to the remote repository of Git for being unable to detect the Git status'; 
    end

end

%% 在matlab命令行中输出执行本函数的结果
disp(statString);
end