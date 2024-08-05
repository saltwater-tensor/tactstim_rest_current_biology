function fh = detectOutputSuppression_tester()
% DETECTOUTPUTSUPPRESSION_TESTER tests the functionality of detectOutputSuppression().
%
%   DETECTOUTPUTSUPPRESSION_TESTER() runs more than 80 tests to confirm the expected
%   outputs and expected errors produced by detectOutputSuppression() for a wide variety
%   of syntaxes. If functioning properly, all intentionally created errors will be
%   bypassed, their error messages will appear on the command window, and a confirmation
%   dialog will appear at the end.
%
%   genericFunction = DETECTOUTPUTSUPPRESSION_TESTER() returns a function handle that
%   should be named "genericFunction" and can be used to run many of the tests from
%   the command window.
%       Example (from the command window):
%           genericFunction = DETECTOUTPUTSUPPRESSION_TESTER();
%           [a,~] = genericFunction([0,1],{'a','~'})
%
%   MATLAB release requirements for these tests.
%   * r2016b and later
%   * A few tests can be commented out so the tester function is fully supported in
%     r2014a and later.
%
%   HOW THIS WORKS
%   The local function genericFunction() is called with various output syntaxes.
%   Within that function, detectOutputSuppression() is called and the results
%   are compared to the expected results which are supplied as inputs to
%   genericFunction().
%
%   Testing expected outputs
%     Example: [a,~] = genericFunction([0,1],{'a','~'})
%   The first input is the expected results of the 1st output to detectOutputSuppression.
%   The second input is the expected results of the 4th output to detectOutputSuppression,
%   components.outNames. Together, this tests that tilde suppression was correctly identified
%   and the outputs were parsed correctly. If outputs do not match expectation, the tester
%   funtion will end with an error message. To stress test the function under many syntaxes,
%   you will see some unorthodox, but accepted MATLAB syntaxes that may require aspirin.
%
%   Testing for expected errors
%     Example: [q{:}, ~] = genericFunction([],{},'DETECTOUTSUP:listExp');
%   Some tests are setup to intentionally cause an error. The 3rd input is the expected
%   error ID. If the expected error is receied, the tester function will continue.  If
%   the error is not received or if there is an error when none was expected, the tester
%   function will end with an error message.
%
%   Much effort was put toward generalizability but that, of course, does not rule out
%   untested syntaxes that may cause problems.  Please report such cases to the developer.
%
% Source: <a href = "https://www.mathworks.com/matlabcentral/fileexchange/79218-detectoutputsuppression">detectOutputSuppression_tester</a>
% Author: <a href = "https://www.mathworks.com/matlabcentral/profile/authors/3753776-adam-danz">Adam Danz</a>
% Copyright (c) 2020  All rights reserved

% Released with vs 2.0 of detectOutputSuppression.

%% Initialize test conditions
if nargout > 0
    fh = @genericFunction;
    return
else
    % Set up tests
    testSetup(true);
end

%% Simple outputs in square brackets with and without suppression
% These syntaxes are all supported by detectOutputSuprression() and should not generate
% any errors from withing detectOutputSuprression.
[a] = genericFunction(0,{'a'});  %#ok<*NASGU>
[~] = genericFunction(1,{'~'});
[a,~] = genericFunction([0,1],{'a','~'}); %#ok<*ASGLU>
[~,b] = genericFunction([1,0],{'~','b'});
[a,~,c] = genericFunction([0,1,0],{'a','~','c'});
[a,b,c] = genericFunction([0,0,0],{'a','b','c'});
[a,~,c,~,e] = genericFunction([0,1,0,1,0],{'a','~','c','~','e'});
[~,b,~,d,~,f] = genericFunction([1,0,1,0,1,0],{'~','b','~','d','~','f'});
[~,b,c,d,e,f,g] = genericFunction([1,0,0,0,0,0,0],{'~','b','c','d','e','f','g'});
[a,b,c,d,e,f,~] = genericFunction([0,0,0,0,0,0,1],{'a','b','c','d','e','f','~'});
[a,~,~,~,~,~,~] = genericFunction([0,1,1,1,1,1,1],{'a','~','~','~','~','~','~'});
[~,~,~,~,~,~,g] = genericFunction([1,1,1,1,1,1,0],{'~','~','~','~','~','~','g'});
[~,~,~,~,~,~,~] = genericFunction([1,1,1,1,1,1,1],{'~','~','~','~','~','~','~'});
[a,b,c,d,e,f,g] = genericFunction([0,0,0,0,0,0,0],{'a','b','c','d','e','f','g'});

% White space should not interfere
[a,~]=genericFunction([0,1],{'a','~'});
[a,~]=   genericFunction([0,1],{'a','~'});
[a,~]    =genericFunction([0,1],{'a','~'});
[a,~] = genericFunction    ([0,1],{'a','~'});
[   a, ~  ] = genericFunction([0,1],{'a','~'});
[~,  b] = genericFunction([1,0],{'~','b'});
[a , ~ , ~ , d , ~ , ~ , g ] = genericFunction([0,1,1,0,1,1,0],{'a','~','~','d','~','~','g'});
    [ ~    ,   b   ]    =     genericFunction  (  [1,0], {'~','b'});  % this line should begin with leading whitespace

%% Outputs not in square brackets
% When using output suppressions, the tilde must be in square brackets.
% These syntaxes are all supported by detectOutputSuprression() and should not
% generate any warnings or errors from withing detectOutputSuprression. However,
% outputs not in brackets will not appear in 'components' output.
a = genericFunction(0,{});
a(2) = genericFunction(0,{});

T = array2table(rand(3,2),'VariableNames', {'A','B'});
T.A(1) = genericFunction(0,{});
c = cell(2,2);
c{2} = genericFunction(0,{});

% No outputs: outputs to detectOutputSuppression() will be empty.
genericFunction({[],{},'',[]}) % This stores inputs for next few lines.
genericFunction();
genericFunction;
genericFunction
fprintf('     ^ the output above was expected.\n')
genericFunction([],{});

% When a function is used as part of an expression, such as an if-statement,
% the 1st output is implicitly assigned and nargout==1. In this case, the first
% output of detectOutputSuppression() will be a scalar logical (0). 
% https://www.mathworks.com/help/matlab/ref/nargout.html#bvj411r-1
genericFunction({0,{},'',[]}) % This stores inputs for the next few lines
genericFunction + 4;
h = genericFunction + 4;
genericFunction()==false; %#ok<EQEFF>
x = genericFunction()==false;
sum(genericFunction());
if isempty(genericFunction()), error('Oops!'), end

% The comment above also applies to these cases when function output
% is implicitly provided as an input when wrapped in another function.
genericFunction({0,{},'',[]}) % This stores inputs for next 2 lines
str= sprintf('First output is %d.', genericFunction());
assert(~isempty(genericFunction()))
eval('x=genericFunction();') % For testing purposes only; don't use eval()!
% Though without outputs within eval(), then no outputs are produced
genericFunction({[],{},'',[]}) % This stores inputs for the next line
eval('genericFunction()') % For testing purposes only; don't use eval()!
fprintf('     ^ the output above was expected.\n')

%% Unrelated content that precedes or follows the call to the inquiry function within the same line
% The parser is able to locate the outputs even when irrelevant material is within the same line.
% These lines are intentionally very messy and examples of really bad sytle.
% These syntaxes are all supported by detectOutputSuprression() and should not
% generate any warnings or errors from withing detectOutputSuprression.
true; [a, b, ~] = genericFunction([0,0,1],{'a','b','~'});
a = rand(); [a, ~, c] = genericFunction([0,1,0],{'a','~','c'});
'    '; [~, b, c] = genericFunction([1,0,0],{'~','b','c'});
disp(' '), [~] = genericFunction(1,{'~'});
disp(' ');, [a] = genericFunction(0,{'a'}); %#ok<*NOCOMMA>
[1,[2:5],[3],[[4,5]]]; [a, b, ~] = genericFunction([0,0,1],{'a','b','~'});   %#ok<*NBRAK>
[a,b] = sort(1:10); [a,b,~,d] = genericFunction([0,0,1,0],{'a','b','~','d'});
[1:5]; [a, b, ~] = genericFunction([0,0,1],{'a','b','~'});
a = 1:5; a([2:3])=0; [~,b] = genericFunction([1,0],{'~','b'}); [a,b] = sort(a);

y = 1:10;
y([1,3,5]); [y,y,~,z] = genericFunction([0,0,1,0],{'y','y','~','z'});  %#ok<*NOEFF>

a = ':'; b = ')'; c = ' ';
[a,b,c], genericFunction([],{}); %#ok<NOPRT,*NOPTS>
fprintf('     ^ the output above was expected.\n')
[a,b,c]; genericFunction([],{}); %#ok<*VUNUS>
[a,b,c]; a = genericFunction(0,{});

%% Variables or functions with very similar names as the inquiry function
% Function or variables with names similar to the inquiry function do not interfere.
% These syntaxes are all supported by detectOutputSuprression() and should not
% generate any warnings or errors from withing detectOutputSuprression.
%
% Variable names that exactly match & shadow the inquiry function would cause a Matlab error.
% Examples:
%   [~, genericFunction] = genericFunction(___);   Error:  'MATLAB:UndefinedFunction'
%   genericFunction = 0; [a]=genericFunction(1);   "genericFunction" is no longer a function.

genericFunction2 = @sort;
[p,q] = genericFunction2(rand(1,5)); [~,~,c,~] = genericFunction([1,1,0,1],{'~','~','c','~'});
[p,q] = genericFunction2(rand(1,5)); [~,~,c,~] = genericFunction ([1,1,0,1],{'~','~','c','~'});

genericFunction({[1,1,0,1],{'~','~','c','~'},'',[]}) % This stores inputs for next few lines
[p,q] = genericFunction2(rand(1,5)); [~,~,c,~] = genericFunction %#ok<NOPRT>
[p,q] = genericFunction2(rand(1,5)); [~,~,c,~] = genericFunction, %#ok<NOPRT>
[p,q] = genericFunction2(rand(1,5)); [~,~,c,~] = genericFunction;
[p,q] = genericFunction2(rand(1,5)); [~,~,c,~] = genericFunction() %#ok<NOPRT>
fprintf('     ^ the 3 outputs above were expected.\n')

genericfunction = num2cell(1:5);
[p,q] = genericfunction{2:3}; [~,~,c,~] = genericFunction([1,1,0,1],{'~','~','c','~'});
[p,q] = genericfunction{[1,5]}; [~,~,c,~] = genericFunction ([1,1,0,1],{'~','~','c','~'});
[p] = genericfunction(1); [~,~,c,~] = genericFunction ([1,1,0,1],{'~','~','c','~'});

%% Other supported syntaxes
% These syntaxes are all supported by detectOutputSuppression() and should not
% generate any warnings or errors from withing detectOutputSuprression.

% More complex output var names
[kcnt32, a_b_c, ~, ~ ,~, ~, A_B_C] = genericFunction([0,0,1,1,1,1,0],{'kcnt32','a_b_c','~','~','~','~','A_B_C'});

% Tables (requires Matlab r2016b or later)
T = array2table(rand(3,2),'VariableNames', {'A','B'});
[T.A(1), ~, T.C(1:3)] = genericFunction([0,1,0],{'T.A','~','T.C'});
[~, T{:,1}, T{1,:}, ~, T{3,2}] = genericFunction([1,0,0,1,0],{'~','T','T','~','T'});

%% Supported indexing within outputs
% detectOutputSuppression() supports scalar output indexing that does not involve comma 
% separated list assignment.
% https://www.mathworks.com/help/matlab/matlab_prog/comma-separated-lists.html.
% These syntaxes are all supported by detectOutputSuprression() and should not
% generate any warnings or errors from withing detectOutputSuprression.

% Indexing with parentheses
[~, b(2), c(3), ~, ~, ~, g(6)] = genericFunction([1,0,0,1,1,1,0],{'~','b','c','~','~','~','g'});
[a(1,2), ~, ~, d(3,1), ~, f(2,1,6), ~] = genericFunction([0,1,1,0,1,0,1],{'a','~','~','d','~','f','~'});
[ a( 1 : 2) , ~, ~, d(3,(1)), ~, f(2,(1:(2))), ~] = genericFunction([0,1,1,0,1,0,1],{'a','~','~','d','~','f','~'});
[a(2:3),b(100:200),~,~,~,f((2*(4-3)):round(5/3)+(8/2))] = genericFunction([0,0,1,1,1,0],{'a','b','~','~','~','f'});

% Scalar indexing with curly brackets
clear a b c d e f g
[~, b{2}, c{3}, ~, ~, ~, g{6}] = genericFunction([1,0,0,1,1,1,0],{'~','b','c','~','~','~','g'});
[a{1,2}, ~, ~, d{3,1}, ~, f{2,1,6}, ~] = genericFunction([0,1,1,0,1,0,1],{'a','~','~','d','~','f','~'});
[~,~,~,~,~,f{2}{1:1}{4}] = genericFunction([1,1,1,1,1,0],{'~','~','~','~','~','f'});
b = cell(1,3);
[~,b{2},~] = genericFunction([1,0,1],{'~','b','~'}); %#ok<*NASGU>

q = cell(1);
[~, q{:}] = genericFunction([1,0],{'~','q'});

q = cell(1,3);
[q{1}, q{2}, ~, q{3}] = genericFunction([0,0,1,0],{'q','q','~','q'});

% Scalar indexing with structures
clear a b c d e f g
[a.a, ~, a.c, ~, a.d.x] = genericFunction([0,1,0,1,0],{'a.a','~','a.c','~','a.d.x'});
[a.a(2), ~, a.c(3,2,4), ~, a.d.x(9:12)] = genericFunction([0,1,0,1,0],{'a.a','~','a.c','~','a.d.x'});

a.a = cell(1,3);
[a.a{1}, ~, a.a{2}, ~, a.a{3}] = genericFunction([0,1,0,1,0],{'a.a','~','a.a','~','a.a'});

S = struct('x',{1,2,cell(1,3),4,5,6});
[S(1).x, S(2).x, S(3).x{2}, ~, ~, S(6).x(99)] = genericFunction([0,0,0,1,1,0],{'S.x','S.x','S.x','~','~','S.x'}); %#ok<*STRNU>

%% Indexing not-supported within outputsOutput (comma separated list assignment)
% Comma separated list assignment is not supported by detectOutputSuppression(). It is detected by
% the discrepancy between the nargout input and the number of ouputs identified by the parser and
% will throw an error 'DETECTOUTSUP:listExp'.  
% https://www.mathworks.com/help/matlab/matlab_prog/comma-separated-lists.html.

% Cell indexing with comma separated list expansion
q = cell(1,3);
[q{:}, ~] = genericFunction([],{},'DETECTOUTSUP:listExp');
[q{1}, q{2:3}, ~] = genericFunction([],{},'DETECTOUTSUP:listExp');

c = 1:3;
a = cell(1,4);
[a{c}] = genericFunction([],{},'DETECTOUTSUP:listExp');

S = struct('x', {1,2,cell(1,2)});
[a, S(2).x, ~, ~, S(3).x{:}, ~] = genericFunction([],{},'DETECTOUTSUP:listExp');

% Structure indexing with comma separated list expansion
S = struct('x',{1,2,3});
[S.x] = genericFunction([],{},'DETECTOUTSUP:listExp');
[~, S(1:2).x] = genericFunction([],{},'DETECTOUTSUP:listExp');
[~, b, S.x] = genericFunction([],{},'DETECTOUTSUP:listExp');
[~, b, S.x] = genericFunction([],{},'DETECTOUTSUP:listExp');

q = cell(1,3);
[q{1}, ~, S.x] = genericFunction([],{},'DETECTOUTSUP:listExp');

%% Expected errors when the nout inputs does not match the actual nargout
% Matlab's nargout will list the correct number of outputs which is why it should be used
% instead of specifying the estimated number of outputs directly in the first input to
% detectOutputSuppression().  If there is a discrepancy between the 'nout' input and the
% number of detected outputs, it will throw an error.

% If the discrepancy was caused by comma separated lists in cell or structure arrays, you'll
% receive the "DETECTOUTSUP:listExp" error ID (demonstrated above).

% When outputs are not separated by commas, the parser will fail and will return "DETECTOUTSUP:noutChk".
[a b] = genericFunction([],{}, 'DETECTOUTSUP:noutChk'); %#ok<NCOMMA>
[a b c ~] = genericFunction([],{}, 'DETECTOUTSUP:noutChk'); %#ok<NCOMMA>
[a, b, c d] = genericFunction([],{}, 'DETECTOUTSUP:noutChk'); %#ok<NCOMMA>
[a, b, c ~] = genericFunction([],{}, 'DETECTOUTSUP:noutChk'); %#ok<NCOMMA>

% If the wrong number of expected outputs was provided in 'nout' (ie, user didn't use nargout),
% the discrepancy will be detected and will throw the "DETECTOUTSUP:noutChk" error ID.
% The 5th input to genericFunction() overrides the nargout input with an incorrect value.
[a,~] = genericFunction([0,1],{'a','~'},'DETECTOUTSUP:noutChk',3);   % nargout==2
[~] = genericFunction(0,{},'DETECTOUTSUP:noutChk',0);                % nargout==1

% If outputs weren't detected and 'nout' does not match nargout, you'll receive a
% 'DETECTOUTSUP:parser' error ID.
genericFunction(0,{},'DETECTOUTSUP:parser',2);                       % nargout==0

%% Expected errors when the inquiry function appears more than once on a line.
% When the inquiry function appears more than once on a line, the parser cannot know which
% line invoked detectOutputSuppression(), even if the function name is commented-out.
% Use an escape character '\' to ignore the function name in a comment (example below).

% Two inquiry functions called on the same line generate the 'multicall' error ID.
[a] = genericFunction([],{},'DETECTOUTSUP:multicall'); [2,3];  [a] = genericFunction([],{},'DETECTOUTSUP:multicall');
[a] = genericFunction([],{},'DETECTOUTSUP:multicall'); [a] = genericFunction([],{},'DETECTOUTSUP:multicall');
[~] = genericFunction([],{},'DETECTOUTSUP:multicall'); genericFunction([],{},'DETECTOUTSUP:multicall');

% This line is OK since the 2nd appearance of the inquiry function is ignored due to the escape character.
[a,~] = genericFunction([0,1],{'a','~'}); % [~,b]=\genericFunction(); OK, no error

%% Expected errors when the inquiry function is wrapped in an anonymous function.
% An inquiry function called from an anonymous function generates the 'noinqfcn' error ID.
% This is because ther parser is expecting an exact match of the inquire function name. 
f = @genericFunction;
j = f([],{},'DETECTOUTSUP:noinqfcn');
[j] = f([],{},'DETECTOUTSUP:noinqfcn');
[~]= f([],{},'DETECTOUTSUP:noinqfcn');
[j,k,l] = f([],{},'DETECTOUTSUP:noinqfcn');
[~,k] = f([],{},'DETECTOUTSUP:noinqfcn');

%% Expected errors when the caller function line is split before the inputs.
% The parser is looking for [outputs]=func and that cannot be split into multiple lines.
% This section will not work with r2014a but works with r2016b.
[a,~,c,...
    ~,e,f] = genericFunction([],{},'DETECTOUTSUP:parser');

[a,b,~,d] = ...
    genericFunction([],{},'DETECTOUTSUP:parser');

% This one below is OK because when nargout indicates 1 output and the parser fails to
% detect square brackets in the line that invoked the inquiry function, it assumes
% the output is not a tilde because tildes must be with square brackets.
a = ...
    genericFunction([0],{});

% This syntax, however, is problematic.  It will *not* generate an error in detectOutputSuppression()
% but will return incorrect results due to the reason explained above!
% % [~] = ...
% %     genericFunction([1],{});

% These are OK since the line split happens after the pattern used by the parser. But these
% tests would cause an error when run from the command winow.
[a,~,~,d] = genericFunction(...
    [0,1,1,0],{'a','~','~','d'});
[a,~,~,d] = genericFunction...
    ([0,1,1,0],{'a','~','~','d'});
[a,~,~,d] = genericFunction([0,1,1,0],...
    {'a','~','~','d'});

%% Expected errors when the input does not match the actual nargout count
% The 'nout' input is expected to be generated by MATLAB's nargout function and is required to
% be a scalar integer greater or equal to 0.  These inputs should cause various Matlab errors.

[~,b] = genericFunction([],{},'MATLAB:narginchk:notEnoughInputs',-1);               % no inputs
[~,b] = genericFunction([],{},'MATLAB:detectOutputSuppression:expectedScalar',-2);  % nout=[];
[~,b] = genericFunction([],{},'MATLAB:detectOutputSuppression:notGreaterEqual',-9); % nout < 0
[~,b] = genericFunction([],{},'MATLAB:detectOutputSuppression:expectedInteger',NaN);% nout not an integer
[~,b] = genericFunction([],{},'MATLAB:detectOutputSuppression:expectedScalar',1:3); % nout not an scalar
[~,b] = genericFunction([],{},'MATLAB:detectOutputSuppression:invalidType','2');    % nout not numeric
[~,b] = genericFunction([],{},'MATLAB:detectOutputSuppression:invalidType',true);   % nout not numeric

%% Confirm success
% All tests have passed.  Generage confirmation dialog. 
testSetup(false)
return

%% Inquiry function
% This is the function inquiring about suppressed outputs.
function [a,b,c,d,e,f,g] = genericFunction(expectedTilde,expectedOutNames,expectedMEID,nout)
% REQUIRED
% * expectedTilde is the expected "isTilde" output (1xn 1s|0s vector; n=nargout); Alternatively it can be a
%   structure that stores the inputs for the next time this func is called. This gives us the ability to test
%   detectOutputSuppression when this func is called without inputs.
% * expectedOutNames is the expected results in components.outNames (1xn cell array)
% OPTIONAL
% * expectedMEID: a string identifying the expected ID for expected errors.
%   When provided, if the error is as expected the code will continue. If
%   the error does not occur or differs from the ID, an error is thrown. If
%   empty or missing, an error is not expected.
% * nout: a scalar value to use in place of nargout as the input to detectOutputSuppression() to cause an
%   error. Default = [], nargout will be used.

% Input validation
persistent nextCallInputs testCounter
if nargin == 0 && ~isempty(nextCallInputs)
    expectedTilde = nextCallInputs{1};
    expectedOutNames = nextCallInputs{2};
    expectedMEID = nextCallInputs{3};
    nout = nextCallInputs{4};
end
if nargin==1 && iscell(expectedTilde)
    nextCallInputs = expectedTilde; % store inputs for next call.
    return
elseif nargin==1 && islogical(expectedTilde) && expectedTilde
    testCounter = 0; % reset test counter
    return
elseif nargin==1 && islogical(expectedTilde) && ~expectedTilde
    a = testCounter; % return test counter
    return
end
if nargin>0 && nargin<3
    expectedMEID = '';
end
if nargin>0 && nargin<4 || isempty(nout)
    % Use to intentionally cause nout error
    nout = nargout;
elseif nout==-2
    nout = []; % to test for missing input.
end
% Test functionality
try
    if nout == -1 % test for no-inputs
        [isTilde, callerTxt, caller, components] = detectOutputSuppression();
    else % all other tests
        [isTilde, callerTxt, caller, components] = detectOutputSuppression(nout);
    end
    % Define outputs
    a=1; b=2; c=3; d=4; e=5; f=6; g=7;
catch ME
    % Define outputs
    a=NaN; b=NaN; c=NaN; d=NaN; e=NaN; f=NaN; g=NaN;
    % Check error
    if isempty(expectedMEID)
        % Unexpected error
        fprintf(2,'\n%s Unexpected error. Error ID: %s\n',char(8999),ME.identifier);
        rethrow(ME)
    elseif ~strcmp(ME.identifier, expectedMEID)
        % Expected error does not match actual error
        fprintf(2,'\n%s Incorrect expected error:\n  Expected error ID: %s\n  Actual error ID: %s\n\n',char(8999),expectedMEID, ME.identifier);
        rethrow(ME)
    else
        % Expected error occured; Show what it would look like in the command window
        % and then leave this local function.
        stackIdx = strcmp({ME.stack.name},mfilename);
        fprintf('%s Test passed for expected error ID ',char(9745))
        fprintf(2,'%s:\n\n%s\n\n',ME.identifier,ME.message);
        if any(stackIdx)
            fprintf('Intentionally caused by %s line %d.\n',ME.stack(stackIdx).name, ME.stack(stackIdx).line)
        else
            fprintf('Intentionally caused.\n') %When called from command window.
        end
        disp([char(95*ones(1,50)),char(10)]) %#ok<CHARTEN> newline() won't work for matlab<16b
        return
    end
end
if ~isempty(expectedMEID)
    error('%s Exected an error with ID "%s" but no error occurred.',char(8999),expectedMEID)
end
% Compare expected-actual results: isTilde
if ~isequal(isTilde,expectedTilde)
    error(['%s Unexpected error in line %d: %s.\n  Expected output suppression index: [%s]\n'...
        '  Computed output suppression index: [%s].'], char(8999),caller.line, callerTxt, ...
        strtrim(regexprep(num2str(expectedTilde),' +',' ')), ...
        strtrim(regexprep(num2str(isTilde),' +',' ')));
end
% Compare expected-actual results: components.outNames
if ~isequal(components.outNames,expectedOutNames)
    error(['%s Unexpected error in line %d: %s.\n  Expected component.outNames: [%s]\n'...
        '  Computed ocomponent.outNames: [%s].'], char(8999),caller.line, callerTxt, ...
        strjoin(expectedOutNames), strjoin(components.outNames));
end
% Increment test counter
testCounter = testCounter+1;

%% Set up tests and indicate success.
function testSetup(TF)
% testSetup(true)
%   Clears command window and documents start of tests. Resets test counter.
% testSetup(false)
%   Displays msgbox indicating success.
headerLine = char(42*ones(1,10));
if TF
    % Initialize command window before tests begin.
    clc(); pause(0.1)
    headerLine = char(42*ones(1,10));
    fprintf('%s Starting %s %s\n',headerLine,mfilename,headerLine)
    % Reset test counter
    genericFunction(true);
else
    % Tests done without unexpected errors; indicate success.
    nTests = genericFunction(false);
    commandwindow()
    msgbox(sprintf(['All %d detectOutputSuppression tests passed. Expected errors are shown in ',...
        'the command window.'],nTests), 'Success')
    fprintf('%s %s complete without unexpected errors. %s\n',headerLine,mfilename,headerLine)
end
