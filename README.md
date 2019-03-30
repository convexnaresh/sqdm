# sqdm

Sqdm is a map data structure for secure point location problem.
It is a structure that secures a map.

##for unicode objects, use No you cannot. Try hashlib.sha1(u'\xcc\x88u') and see. – Sassa NF Mar 13 '17 at 10:50
##@SassaNF that's because the u'...' causes the string to be a unicode object instead of a str object.
## To turn a unicode object into a str, you'll need to
# encode it. e.g. hashlib.sha1(u'\xcc\x88u'.encode("utf-8"))

