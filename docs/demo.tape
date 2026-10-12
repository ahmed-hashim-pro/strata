# Terminal demo for the README. Regenerate with:  make demo
# Records the real CLI against examples/site.csv; nothing is staged.

Output docs/demo.gif

Set Shell "bash"
Set FontSize 15
Set Width 1060
Set Height 520
Set Padding 18
Set Theme "GitHub Dark"
Set TypingSpeed 35ms

Hide
Type "export PATH=/opt/strata:$PATH PS1='$ ' && clear"
Enter
Show

Type "strata query examples/site.csv 490 440    # next to hole BH16"
Enter
Sleep 3s

Type "strata query examples/site.csv 600 480    # 120 m from the nearest hole"
Enter
Sleep 4s

Type "echo exit $?"
Enter
Sleep 1.5s

Type "clear"
Enter
Sleep 300ms

Type "strata bench --sites 50    # 50 synthetic sites, checked against the true geology"
Enter
Sleep 7s

Type "clear"
Enter
Sleep 300ms

Type "strata validate tests/data/typo.csv"
Enter
Sleep 3s
