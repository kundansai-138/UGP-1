for i in $(seq 1 200); do
  ping -c1 -W1 csews$i >/dev/null 2>&1 && echo "csews$i up"
done | tee up_nodes.txt
