function statString = configureKits(statNum, msg)
% ��;��
% 1�����������ݻ�ǰ��ʹ�øú������л������ã��ر����д��ڣ�������������ݣ�Ȼ����Ӹ���Toolkit��·����
% 2���������ݻ��㷨����ʱ��ͨ�����ñ�������ʵ�ֶ������ݻ����߰���ʵʱ�汾������Զ�̴����ύ
% Ҫ���û�����ʹ�ñ������ṩ�İ汾���ƹ��ܣ���Ҫ�����ڱ�������Git�Ŀͻ��ˣ����ӣ�http://git-scm.com/download/win��
% ������ɺ�ȷ��Git��ִ���ļ���·���Ѿ�����ӵ���ϵͳ·����Path���У�һ����԰�װʱ�ǻ�Ĭ����ӵ�
% �˺�����ȷ������Git�Ļ����ϣ����ܵ��ñ�����
% 
% �������
% statNum��int���ͱ���������������
% ��statNumΪ0����ʾ���øú�����������ݻ��Ļ�����������
% ��statNumΪ1����ʾ���øú�����Git��Զ�̴���ֿ⣨��Git�������ˣ�
% ��statNumΪ2����ʾ���øú��������е��޸�ȫ����ӵ�Git�ı����ݴ���
% ��statNumΪ3����ʾ���øú��������е��޸�ȫ���ύ��Git�ı��ش���ֿ�
% ��statNumΪ4����ʾ���øú��������е��޸�ȫ���ύ��Git��Զ�̴���ֿ�
% msg���ַ������ͱ�������statNumΪ3ʱ�����ʹ��msg����Ե�ǰ�����޸ĵ�ע�ͣ����û�û������msg���������ύĬ��ע��
% �������
% statString���ַ������ͱ�����configureKits������ִ�н��
%
% ���ޣ�����Ŀǰ��������Windows����ϵͳ��Linux�汾��δ����

%% ������������
try
    statusNum = statNum;
catch  % ����û�û������״̬��
    statusNum = 0; % Ĭ��Ϊ0
end

%% ��statusNumΪ0���ú����������û����境
if statusNum == 0
    close all
    clc;
    format long g
    addpath('./CCSToolkit'); % ���CCS Toolkit���߰�
    addpath('./EvolveToolkit'); % ���Evolve Toolkit���߰�
    addpath('./CCSDebugKit'); % ���CCS�ݻ���Debug���߰�
    statString = 'Environment configuration has been completed';

%% ��statusNumΪ1����Git��Զ�ֿ̲����ش��뵽�������ڸ���
elseif statusNum == 1
    [statCode, resultString] = system('git pull');
    if statCode == 0
         % ������صĺ���ִ�н���ַ����а���'Already up-to-date',˵�����ش���İ汾�Ѿ���Git�������˵����°汾
         % ���߱�Git�������˵İ汾��Ҫ�£�Ҳ�����ڱ����޸Ļ�û�ύ���������˵Ĵ���
        if ~isempty(strfind(resultString, 'Already up-to-date'))
            statString = 'The code in this computer is already the newest version';
        % ������صĺ���ִ�н���ַ����в�����'Already up-to-date',��˵��������Ҫ���£�����ͨ��git
        % pull������µ������°汾
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

%% ��statusNumΪ2���ύ�޸ĵ�Git�ı����ݴ���    
elseif statusNum == 2
    [~, resultString] = system('git status');
    checkStr1 = 'Changes not staged for commit:';
    checkStr2 = 'Changes to be committed:';
    checkStr3 = 'nothing to commit, working directory clean';
    checkStr4 = 'You have unmerged paths.';
    checkStr5 = 'nothing added to commit but untracked files present';
    % ������صĺ���ִ�н���ַ����а���checkStr1������,˵�����޸���δ�ύ���ݴ�������ʱ���ܵ���git add -A
    if ~isempty(strfind(resultString, checkStr1)) 
        [statCode, ~] = system('git add -A');
        if statCode == 0 % ���״̬��Ϊ0��������ȷ�ύ
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
        if statCode == 0 % ���״̬��Ϊ0��������ȷ�ύ
            statString = 'New Files detected and added to the local stage area';
        else
            statString = 'Failed to add the new files to the local stage area';
        end
    else
        statString = 'Failed to add the changes to the local stage area for being unable to detect the Git status'; 
    end
    
%% ��statusNumΪ3���ύ�޸ĵ�Git�ı��زֿ�
elseif statusNum == 3
    try
        msgString = msg;
    catch  % ����û�û�����뽫�����ύ�����زֿ�ʱ��ע��
        msgString = 'new modifications'; % Ĭ��ע��Ϊ��new modifications��
    end
    cmdString = ['git commit -m "' msgString '"'];
    [statCode, resultString] = system(cmdString);
    if statCode == 0 % ���״̬��Ϊ0��������ȷ�ύ
        statString = 'The changes of code have been committed to the local repository of Git';
    else % ���״̬�벻Ϊ0
        checkStr1 = 'nothing to commit, working directory clean';
        checkStr2 = 'Changes not staged for commit:';
        checkStr3 = 'fatal: Exiting because of an unresolved conflict.';
        % ������صĺ���ִ�н���ַ����а���checkStr1������,˵�������������ύ�������ύ����
        if ~isempty(strfind(resultString, checkStr1)) 
            statString = 'Nothing new to commit to the local repository of Git, working directory is clean';
        % ������صĺ���ִ�н���ַ����а���checkStr2������,˵���û���δ����git add�����ύ�޸ĵ��ݴ�������ʱ��Ҫ�ȵ���configureKits(2)
        elseif ~isempty(strfind(resultString, checkStr2)) 
            statString = 'Failed to commit the changes, please use "configureKits(2)" to add the changes to the stage area first';
        elseif ~isempty(strfind(resultString, checkStr3)) 
            statString = 'Failed to commit the changes, please use "configureKits(2)" to fix the conflicts of merge first';
        else % ���������checkStr1��checkStr2�����ݣ�˵���ύ����
            statString = 'Failed to commit the changes to the local repository of Git for unknown reason';
        end
    end

%% ��statusNumΪ4���ύ�޸ĵ�Git��Զ�ֿ̲�
elseif statusNum == 4
    % ����Git���������û���push֮ǰ�Ƿ��Ѿ��ڱ���add��commit��������ύ���뵽Զ�ֿ̲�ǰ����Ҫ���ȼ���ǰ���޸��Ƿ��Ѿ��ύ
    [~, resultString] = system('git status');
    checkStr1 = 'Changes not staged for commit:';
    checkStr2 = 'Changes to be committed:';
    checkStr3 = 'nothing to commit, working directory clean';
    % ������صĺ���ִ�н���ַ����а���checkStr1������,˵�����������޸���δ��ӵ��ݴ���
    if ~isempty(strfind(resultString, checkStr1))
        statString = 'Failed to commit the changes, please use "configureKits(2)" to add the changes to the stage area first';
    % ������صĺ���ִ�н���ַ����а���checkStr2������,˵�����е��޸ľ�����ӵ��ݴ�����������δ�ύ�����زֿ�
    elseif ~isempty(strfind(resultString, checkStr2))
        statString = 'Failed to commit the changes, please use "configureKits(3, your notes)" to add the changes to the local repository first';
    % ������صĺ���ִ�н���ַ����а���checkStr3������,˵�����е��޸����ύ�����زֿ�
    elseif ~isempty(strfind(resultString, checkStr3))
    % ͨ����������֤���Ż����git pull���д����Զ���ύ
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

%% ��matlab�����������ִ�б������Ľ��
disp(statString);
end