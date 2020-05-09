item = 'SRR8368689'
with open('/tmp/' + item + '_old.blast') as f:
    matches = f.read()
    for k in range(1, 11):
        print(matches.count('espaceur' + 'old'.capitalize() + str(k) + ','))